# =============================================================================
# TrackMate Analysis — Comprehensive Per-Species 4D Embryo Nuclei Tracking
# =============================================================================
#
# PURPOSE:
#   Step 3 of the production pipeline:
#     Step 0: TrackMate → export ALL tracks (no filters)
#     Step 1: embryo_viewer.py → orient, sphere-fit, ROI, ingression detection
#     Step 2: trackmate_filter_and_validate.R → filter + validate (medaka only)
#   → Step 3: THIS SCRIPT → comprehensive per-species biological analysis
#
# DESIGN:
#   • SEPARATE per-species PDFs (not combined) — each species gets its own set
#   • Cross-species comparison dashboard as a final integrative section
#   • Uses enriched columns from napari viewer: radial/tangential velocity,
#     spherical depth, ingression labels, ROI flags
#   • Computes instantaneous 3D velocity, MSD, turning angles from raw positions
#   • Comprehensive ingression analysis for zebrafish (validated labels)
#   • Margin/ROI density analysis for ingression characterisation
#
# INPUTS (per species, from analysis_output_{species}/):
#   - oriented_spots.csv      38-col enriched spot data
#   - track_summary.csv       20-col per-track aggregates
#   - sphere_params.csv       Fitted sphere centre & radius
#   - roi_bounds.csv          ROI bounding box
#   - ingression_params.csv   Ingression detection parameters
#   - analysis_metadata.json  Viewer session metadata
#   For medaka (filtered): analysis_output/filtered_spots.csv, filtered_track_summary.csv
#
# OUTPUTS (in analysis_output/):
#   Per species ({zf,mk}_):
#     01_overview_qc.pdf        — nuclei count, rate of change, summary
#     02_spatial.pdf            — theta-phi, XY, depth, angular coverage
#     03_track_topology.pdf     — length, duration, birth/death, active
#     04_density.pdf            — 3D voxel density, ROI density
#     05_velocity.pdf           — radial, tangential, instantaneous speed
#     06_flow_fields.pdf        — velocity vector fields (theta-phi, XY)
#     07_physical_properties.pdf — MSD, confinement, displacement, diffusion
#     08_directionality.pdf     — inward fraction, sustained, turning angles
#     09_intensity_qc.pdf       — intensity, SNR, quality
#   Zebrafish only:
#     zf_10_ingression.pdf      — onset timing, spatial, pseudo-trajectories
#     zf_11_margin_density.pdf  — cell density in/near ingression margin
#   Cross-species:
#     12_cross_species.pdf      — normalised comparisons
#   Data:
#     analysis_summary.csv
# =============================================================================

source("renv/activate.R")

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(viridis)
library(purrr)
library(tibble)
library(patchwork)
library(scales)
library(data.table)
library(jsonlite)

# =============================================================================
# PARAMETERS
# =============================================================================

FRAME_INTERVAL_SEC <- 30
FRAME_INTERVAL_MIN <- FRAME_INTERVAL_SEC / 60
SPEED_CONVERSION   <- 60 / FRAME_INTERVAL_SEC   # µm/frame → µm/min

OUTPUT_DIR <- "analysis_output"

SPECIES <- list(
  zebrafish = list(
    name      = "Zebrafish",
    short     = "ZF",
    prefix    = "zf",
    color     = "#2166AC",
    data_dir  = "analysis_output_zebrafish",
    has_ingression = TRUE
  ),
  medaka = list(
    name      = "Medaka",
    short     = "MK",
    prefix    = "mk",
    color     = "#B2182B",
    data_dir  = "analysis_output_medaka",
    filter_dir = "analysis_output",
    has_ingression = FALSE
  )
)

# Voxel density
VOXEL_SIZE_UM <- 50

# Depth bins (µm from surface)
DEPTH_BREAKS <- c(-Inf, 0, 10, 20, 40, 60, Inf)
DEPTH_LABELS <- c("Outside (<0)", "Surface (0-10)", "Shallow (10-20)",
                   "Mid (20-40)", "Deep (40-60)", "Very deep (>60)")

# MSD analysis
MSD_MAX_LAG    <- 30   # max lag in frames
MSD_N_TRACKS   <- 8000 # max tracks to sample for MSD

# Velocity field grid
FLOW_BIN_THETA <- 10   # degrees
FLOW_BIN_PHI   <- 10
FLOW_BIN_XY    <- 30   # µm
FLOW_MIN_N     <- 20   # minimum spots per bin for flow field

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# =============================================================================
# THEME & HELPERS
# =============================================================================

species_colors <- c("Zebrafish" = "#2166AC", "Medaka" = "#B2182B")

theme_pipe <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      plot.subtitle    = element_text(hjust = 0.5, size = base_size - 1, color = "grey50"),
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold")
    )
}

save_pdf <- function(plot, filename, width = 12, height = 8) {
  path <- file.path(OUTPUT_DIR, filename)
  ggsave(path, plot, width = width, height = height, device = cairo_pdf)
  cat(sprintf("  -> Saved %s\n", path))
}

frame_to_min <- function(f) f * FRAME_INTERVAL_MIN
frame_to_hpf <- function(f) f * FRAME_INTERVAL_SEC / 3600

add_depth_bin <- function(df) {
  df %>%
    mutate(depth_bin = cut(SPHERICAL_DEPTH,
                           breaks = DEPTH_BREAKS,
                           labels = DEPTH_LABELS,
                           ordered_result = TRUE))
}

# =============================================================================
# HELPER: compute instantaneous 3D velocity from consecutive spots
# =============================================================================

compute_inst_velocity <- function(spots_df) {
  spots_df %>%
    arrange(TRACK_ID, FRAME) %>%
    group_by(TRACK_ID) %>%
    mutate(
      dx = POSITION_X - lag(POSITION_X),
      dy = POSITION_Y - lag(POSITION_Y),
      dz = POSITION_Z - lag(POSITION_Z),
      dt_frames = FRAME - lag(FRAME),
      displacement_3d = sqrt(dx^2 + dy^2 + dz^2),
      inst_speed_um_min = displacement_3d / (dt_frames * FRAME_INTERVAL_MIN),
      # Component velocities (um/frame)
      vx = dx / dt_frames,
      vy = dy / dt_frames,
      vz = dz / dt_frames,
      # Tangential speed: |v_total|^2 = v_radial^2 + v_tangential^2
      total_speed_frame = displacement_3d / dt_frames,
      tangential_speed = sqrt(pmax(total_speed_frame^2 -
                                     ifelse(is.na(RADIAL_VELOCITY), 0,
                                            RADIAL_VELOCITY)^2, 0)),
      # Turning angle (degrees) between consecutive displacement vectors
      dot_prev = dx * lag(dx) + dy * lag(dy) + dz * lag(dz),
      mag_curr = displacement_3d,
      mag_prev = lag(displacement_3d),
      cos_angle = dot_prev / (mag_curr * mag_prev),
      turning_angle = acos(pmin(pmax(cos_angle, -1), 1)) * 180 / pi,
      # XY direction angle (degrees from +X axis)
      direction_xy = atan2(dy, dx) * 180 / pi
    ) %>%
    ungroup()
}

# =============================================================================
# HELPER: compute MSD for tracks
# =============================================================================

compute_msd_ensemble <- function(spots_df, max_lag = 30, n_sample = 8000) {
  set.seed(42)
  track_ids <- unique(spots_df$TRACK_ID)
  if (length(track_ids) > n_sample) {
    track_ids <- sample(track_ids, n_sample)
  }

  dt <- as.data.table(spots_df)[TRACK_ID %in% track_ids,
                                 .(TRACK_ID, FRAME, POSITION_X, POSITION_Y, POSITION_Z)]
  setkey(dt, TRACK_ID, FRAME)

  msd_list <- list()
  for (lag in 1:max_lag) {
    d2 <- dt[, {
      n <- .N
      if (n > lag) {
        idx1 <- 1:(n - lag)
        idx2 <- (1 + lag):n
        # Only use pairs where frame difference matches lag
        valid <- (FRAME[idx2] - FRAME[idx1]) == lag
        if (sum(valid) > 0) {
          dsq <- (POSITION_X[idx2[valid]] - POSITION_X[idx1[valid]])^2 +
                 (POSITION_Y[idx2[valid]] - POSITION_Y[idx1[valid]])^2 +
                 (POSITION_Z[idx2[valid]] - POSITION_Z[idx1[valid]])^2
          .(msd = mean(dsq), n = length(dsq))
        } else {
          .(msd = numeric(0), n = integer(0))
        }
      } else {
        .(msd = numeric(0), n = integer(0))
      }
    }, by = TRACK_ID]

    if (nrow(d2) > 0 && sum(d2$n) > 0) {
      msd_list[[lag]] <- data.frame(
        lag_frames = lag,
        lag_min    = lag * FRAME_INTERVAL_MIN,
        mean_msd   = weighted.mean(d2$msd, d2$n),
        n_pairs    = sum(d2$n),
        n_tracks   = nrow(d2)
      )
    }
  }
  bind_rows(msd_list)
}

# #############################################################################
#
#                         MAIN EXECUTION
#
# #############################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("  TRACKMATE COMPREHENSIVE ANALYSIS (per-species)\n")
cat(strrep("=", 70), "\n")
cat(sprintf("  Frame interval: %d sec | Speed conversion: x%.0f\n",
            FRAME_INTERVAL_SEC, SPEED_CONVERSION))
cat(sprintf("  Output: %s/\n", OUTPUT_DIR))

# #############################################################################
# SECTION 00: DATA LOADING
# #############################################################################

cat("\n=== SECTION 00: LOADING DATA ===\n")

load_species_data <- function(sp) {
  dir <- sp$data_dir
  if (!dir.exists(dir)) {
    warning(sprintf("Data directory '%s' not found -- skipping %s", dir, sp$name))
    return(NULL)
  }

  cat(sprintf("\n  Loading %s from %s/ ...\n", sp$name, dir))

  # --- Spots (prefer filtered version if filter_dir is set) ---
  filtered_spots <- NULL
  filtered_ts    <- NULL
  if (!is.null(sp$filter_dir)) {
    fp <- file.path(sp$filter_dir, "filtered_spots.csv")
    ft <- file.path(sp$filter_dir, "filtered_track_summary.csv")
    if (file.exists(fp)) filtered_spots <- fp
    if (file.exists(ft)) filtered_ts    <- ft
  }

  if (!is.null(filtered_spots)) {
    spots_path <- filtered_spots
    cat(sprintf("    * Using FILTERED spots: %s\n", spots_path))
  } else {
    spots_path <- file.path(dir, "oriented_spots.csv")
  }

  spots <- fread(spots_path, showProgress = FALSE)
  setDF(spots)
  spots$species <- sp$name

  # Ensure INGRESSING is logical
  if ("INGRESSING" %in% names(spots)) {
    spots$INGRESSING <- as.logical(spots$INGRESSING)
  }
  if ("IN_ROI" %in% names(spots)) {
    spots$IN_ROI <- as.logical(spots$IN_ROI)
  }

  cat(sprintf("    Spots: %s rows, %d columns, frames %d-%d\n",
              format(nrow(spots), big.mark = ","), ncol(spots),
              min(spots$FRAME), max(spots$FRAME)))

  # --- Track summary (prefer filtered version) ---
  if (!is.null(filtered_ts)) {
    ts_path <- filtered_ts
    cat(sprintf("    * Using FILTERED track summary: %s\n", ts_path))
  } else {
    ts_path <- file.path(dir, "track_summary.csv")
  }

  tsummary <- read_csv(ts_path, show_col_types = FALSE)
  tsummary$species <- sp$name
  if ("INGRESSING" %in% names(tsummary)) {
    tsummary$INGRESSING <- as.logical(tsummary$INGRESSING)
  }
  if ("IN_ROI" %in% names(tsummary)) {
    tsummary$IN_ROI <- as.logical(tsummary$IN_ROI)
  }
  cat(sprintf("    Tracks: %s\n", format(nrow(tsummary), big.mark = ",")))

  # --- Auxiliary data ---
  sphere <- read_csv(file.path(dir, "sphere_params.csv"), show_col_types = FALSE)
  sphere_list <- setNames(sphere$value, sphere$parameter)

  roi <- read_csv(file.path(dir, "roi_bounds.csv"), show_col_types = FALSE)
  roi_list <- setNames(roi$value, roi$parameter)

  ingr <- read_csv(file.path(dir, "ingression_params.csv"), show_col_types = FALSE)
  ingr_list <- setNames(as.character(ingr$value), ingr$parameter)

  meta <- fromJSON(file.path(dir, "analysis_metadata.json"))

  list(
    spots      = spots,
    tracks     = tsummary,
    sphere     = sphere_list,
    roi        = roi_list,
    ingression = ingr_list,
    metadata   = meta,
    config     = sp
  )
}

# Load all species
DATA <- list()
for (sp_key in names(SPECIES)) {
  result <- load_species_data(SPECIES[[sp_key]])
  if (!is.null(result)) DATA[[sp_key]] <- result
}

if (length(DATA) == 0) stop("No species data loaded. Check data directories.")

cat(sprintf("\n  Loaded %d species: %s\n",
            length(DATA),
            paste(sapply(DATA, function(d) d$config$name), collapse = ", ")))

# --- Data overview table ---
overview_rows <- lapply(DATA, function(d) {
  sp <- d$spots; ts <- d$tracks
  tibble(
    Species              = d$config$name,
    `Total spots`        = nrow(sp),
    `Total tracks`       = nrow(ts),
    `Frames`             = length(unique(sp$FRAME)),
    `Frame range`        = sprintf("%d - %d", min(sp$FRAME), max(sp$FRAME)),
    `Duration (min)`     = round(max(sp$FRAME) * FRAME_INTERVAL_MIN, 1),
    `Sphere R (um)`      = round(as.numeric(d$sphere["radius"]), 1),
    `Median track length` = median(ts$n_spots),
    `Mean spots/frame`   = round(nrow(sp) / length(unique(sp$FRAME)), 0),
    `Depth range (um)`   = sprintf("%.1f - %.1f",
                                    min(sp$SPHERICAL_DEPTH, na.rm = TRUE),
                                    max(sp$SPHERICAL_DEPTH, na.rm = TRUE)),
    `Ingression validated` = d$config$has_ingression,
    `N ingressing tracks`  = if (d$config$has_ingression) sum(ts$INGRESSING, na.rm = TRUE) else NA_integer_
  )
})
overview_df <- bind_rows(overview_rows)
print(as_tibble(overview_df))
write_csv(overview_df, file.path(OUTPUT_DIR, "analysis_summary.csv"))
cat("  -> Saved analysis_summary.csv\n")

# --- Compute instantaneous velocities for all species ---
cat("\n  Computing instantaneous velocities ...\n")
for (sp_key in names(DATA)) {
  cat(sprintf("    %s ...\n", DATA[[sp_key]]$config$name))
  DATA[[sp_key]]$spots_vel <- compute_inst_velocity(DATA[[sp_key]]$spots)
}

# #############################################################################
#
#          PER-SPECIES ANALYSIS (loop over each species)
#
# #############################################################################

for (sp_key in names(DATA)) {
  d       <- DATA[[sp_key]]
  sp      <- d$config
  prefix  <- sp$prefix
  spots   <- d$spots
  tracks  <- d$tracks
  spots_v <- d$spots_vel
  R_sphere <- as.numeric(d$sphere["radius"])
  sp_col  <- sp$color

  cat(sprintf("\n%s\n  PER-SPECIES ANALYSIS: %s\n%s\n",
              strrep("=", 70), sp$name, strrep("=", 70)))

  # =========================================================================
  # S01: OVERVIEW & QC
  # =========================================================================

  cat(sprintf("\n=== %s S01: OVERVIEW & QC ===\n", sp$short))

  # Nuclei count per frame
  count_df <- spots %>%
    group_by(FRAME) %>%
    summarise(n_nuclei = n(), .groups = "drop") %>%
    mutate(time_min = frame_to_min(FRAME))

  p_count <- ggplot(count_df, aes(x = time_min, y = n_nuclei)) +
    geom_line(color = sp_col, linewidth = 0.5, alpha = 0.7) +
    geom_smooth(method = "loess", span = 0.3, se = TRUE, color = "black",
                linewidth = 0.8) +
    labs(title = sprintf("Detected Nuclei Over Time -- %s", sp$name),
         subtitle = sprintf("Total: %s spots, %s tracks, %d frames",
                            format(nrow(spots), big.mark = ","),
                            format(nrow(tracks), big.mark = ","),
                            length(unique(spots$FRAME))),
         x = "Time (min)", y = "Number of nuclei") +
    theme_pipe()

  # Rate of change
  rate_df <- count_df %>%
    arrange(FRAME) %>%
    mutate(delta = n_nuclei - lag(n_nuclei),
           pct_change = delta / lag(n_nuclei) * 100)

  p_rate <- rate_df %>%
    filter(!is.na(delta)) %>%
    ggplot(aes(x = time_min, y = delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(alpha = 0.3, color = sp_col, linewidth = 0.3) +
    geom_smooth(method = "loess", span = 0.2, se = FALSE,
                color = sp_col, linewidth = 0.8) +
    labs(title = "Frame-to-Frame Change in Nuclei Count",
         subtitle = "Positive = gain (division/detection), Negative = loss (tracking dropout)",
         x = "Time (min)", y = "Delta Nuclei / frame") +
    theme_pipe()

  # Data summary as text annotation plot
  info_text <- sprintf(
    paste0("Species: %s\nSphere R: %.1f um\nTotal spots: %s\nTotal tracks: %s\n",
           "Frames: %d (%s - %s)\nDuration: %.1f min\n",
           "Median track length: %d spots\nMean spots/frame: %.0f\n",
           "Depth range: %.1f - %.1f um%s"),
    sp$name, R_sphere,
    format(nrow(spots), big.mark = ","),
    format(nrow(tracks), big.mark = ","),
    length(unique(spots$FRAME)),
    min(spots$FRAME), max(spots$FRAME),
    max(spots$FRAME) * FRAME_INTERVAL_MIN,
    median(tracks$n_spots),
    nrow(spots) / length(unique(spots$FRAME)),
    min(spots$SPHERICAL_DEPTH, na.rm = TRUE),
    max(spots$SPHERICAL_DEPTH, na.rm = TRUE),
    if (sp$has_ingression) sprintf("\nIngressing tracks: %d (%.1f%%)",
        sum(tracks$INGRESSING, na.rm = TRUE),
        100 * sum(tracks$INGRESSING, na.rm = TRUE) / nrow(tracks)) else ""
  )

  p_info <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = info_text,
             hjust = 0.5, vjust = 0.5, size = 3.5, family = "mono") +
    theme_void() +
    labs(title = "Dataset Summary") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

  p01 <- p_count / p_rate / p_info +
    plot_annotation(
      title = sprintf("%s -- 01 Overview & QC", sp$name),
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

  save_pdf(p01, sprintf("%s_01_overview_qc.pdf", prefix), width = 14, height = 14)


  # =========================================================================
  # S02: SPATIAL DISTRIBUTION
  # =========================================================================

  cat(sprintf("\n=== %s S02: SPATIAL DISTRIBUTION ===\n", sp$short))

  # Theta-phi coverage heatmap
  p_angular <- spots %>%
    filter(!is.na(THETA_DEG), !is.na(PHI_DEG)) %>%
    ggplot(aes(x = PHI_DEG, y = THETA_DEG)) +
    geom_bin2d(bins = c(60, 40)) +
    scale_fill_viridis_c(option = "plasma", trans = "log10", labels = comma) +
    labs(title = "Angular Coverage (Spherical Coordinates)",
         subtitle = "Theta = polar angle from AP, Phi = azimuthal",
         x = "Phi (azimuthal, deg)", y = "Theta (polar, deg)", fill = "Count") +
    theme_pipe()

  # XY projection density
  p_xy <- spots %>%
    ggplot(aes(x = POSITION_X, y = POSITION_Y)) +
    geom_bin2d(bins = 80) +
    scale_fill_viridis_c(option = "plasma", trans = "log10", labels = comma) +
    labs(title = "XY Projection of All Nuclei",
         x = "X (um)", y = "Y (um)", fill = "Count") +
    theme_pipe()

  # Spherical depth distribution
  p_depth_hist <- spots %>%
    ggplot(aes(x = SPHERICAL_DEPTH)) +
    geom_histogram(bins = 60, fill = sp_col, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    labs(title = "Spherical Depth Distribution",
         subtitle = "0 = fitted surface, positive = inward",
         x = "Spherical depth (um)", y = "Count") +
    theme_pipe()

  # Depth x Theta heatmap
  p_depth_theta <- spots %>%
    filter(!is.na(THETA_DEG)) %>%
    ggplot(aes(x = THETA_DEG, y = SPHERICAL_DEPTH)) +
    geom_bin2d(bins = c(50, 50)) +
    scale_fill_viridis_c(option = "inferno", trans = "log10", labels = comma) +
    labs(title = "Depth vs. Polar Angle",
         x = "Theta (polar, deg)", y = "Spherical depth (um)", fill = "Count") +
    theme_pipe()

  # Depth over time heatmap
  p_depth_time <- spots %>%
    mutate(time_min = frame_to_min(FRAME)) %>%
    ggplot(aes(x = time_min, y = SPHERICAL_DEPTH)) +
    geom_bin2d(bins = c(80, 60)) +
    scale_fill_viridis_c(option = "inferno", trans = "log10", labels = comma) +
    labs(title = "Depth x Time Heatmap",
         x = "Time (min)", y = "Spherical depth (um)", fill = "Count") +
    theme_pipe()

  # Track centroids colored by depth
  p_centroids <- tracks %>%
    ggplot(aes(x = mean_x, y = mean_y, color = mean_SPHERICAL_DEPTH)) +
    geom_point(alpha = 0.15, size = 0.3) +
    scale_color_viridis_c(option = "mako", direction = -1) +
    labs(title = "Track Centroid Positions (XY)",
         subtitle = "Coloured by mean spherical depth",
         x = "Mean X (um)", y = "Mean Y (um)", color = "Depth (um)") +
    theme_pipe()

  p02 <- (p_angular | p_xy) / (p_depth_hist | p_depth_theta) / (p_depth_time | p_centroids) +
    plot_annotation(
      title = sprintf("%s -- 02 Spatial Distribution", sp$name),
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

  save_pdf(p02, sprintf("%s_02_spatial.pdf", prefix), width = 16, height = 16)


  # =========================================================================
  # S03: TRACK TOPOLOGY
  # =========================================================================

  cat(sprintf("\n=== %s S03: TRACK TOPOLOGY ===\n", sp$short))

  # Track length distribution
  p_tracklen <- tracks %>%
    ggplot(aes(x = n_spots)) +
    geom_histogram(bins = 80, fill = sp_col, alpha = 0.7) +
    scale_x_log10(labels = comma) +
    labs(title = "Track Length Distribution",
         subtitle = sprintf("Median: %d spots", median(tracks$n_spots)),
         x = "Track length (spots, log)", y = "Count") +
    theme_pipe()

  # Track temporal coverage
  p_coverage <- tracks %>%
    mutate(duration_min = frame_to_min(frame_end - frame_start),
           coverage_pct = n_spots / (frame_end - frame_start + 1) * 100) %>%
    ggplot(aes(x = duration_min, y = coverage_pct)) +
    geom_point(alpha = 0.03, size = 0.3, color = sp_col) +
    geom_density2d(linewidth = 0.4, color = "black") +
    coord_cartesian(ylim = c(0, 105)) +
    labs(title = "Duration vs. Temporal Coverage",
         subtitle = "Coverage = spots / expected frames x 100%",
         x = "Duration (min)", y = "Coverage (%)") +
    theme_pipe()

  # Birth/death timing
  track_timing <- tracks %>%
    select(frame_start, frame_end) %>%
    pivot_longer(everything(), names_to = "event", values_to = "frame") %>%
    mutate(event = ifelse(event == "frame_start", "Track starts", "Track ends"),
           time_min = frame_to_min(frame))

  p_timing <- track_timing %>%
    ggplot(aes(x = time_min, fill = event)) +
    geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c("Track starts" = "#66C2A5", "Track ends" = "#FC8D62")) +
    labs(title = "Track Birth / Death Over Time",
         x = "Time (min)", y = "Count", fill = "") +
    theme_pipe()

  # Active tracks over time
  frames_seq <- seq(min(tracks$frame_start), max(tracks$frame_end))
  active_df <- tibble(
    FRAME = frames_seq,
    n_active = sapply(frames_seq, function(f) {
      sum(tracks$frame_start <= f & tracks$frame_end >= f)
    }),
    time_min = frame_to_min(frames_seq)
  )

  p_active <- ggplot(active_df, aes(x = time_min, y = n_active)) +
    geom_line(color = sp_col, linewidth = 0.6) +
    labs(title = "Active Tracks Over Time",
         x = "Time (min)", y = "Active tracks") +
    theme_pipe()

  p03 <- (p_tracklen | p_coverage) / (p_timing | p_active) +
    plot_annotation(
      title = sprintf("%s -- 03 Track Topology", sp$name),
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

  save_pdf(p03, sprintf("%s_03_track_topology.pdf", prefix), width = 14, height = 10)


  # =========================================================================
  # S04: NUCLEAR DENSITY
  # =========================================================================

  cat(sprintf("\n=== %s S04: NUCLEAR DENSITY ===\n", sp$short))

  # 3D voxel density
  dt <- as.data.table(spots)[, .(FRAME, POSITION_X, POSITION_Y, POSITION_Z, IN_ROI)]
  dt[, `:=`(vx = floor(POSITION_X / VOXEL_SIZE_UM),
            vy = floor(POSITION_Y / VOXEL_SIZE_UM),
            vz = floor(POSITION_Z / VOXEL_SIZE_UM))]

  voxel_counts <- dt[, .(count = .N), by = .(FRAME, vx, vy, vz)]

  dens_stats <- voxel_counts[, .(
    mean_density   = mean(count),
    median_density = as.double(median(count)),
    max_density    = max(count),
    q25            = quantile(count, 0.25),
    q75            = quantile(count, 0.75)
  ), by = FRAME]
  dens_stats[, time_min := frame_to_min(FRAME)]

  p_dens_mean <- ggplot(as.data.frame(dens_stats), aes(x = time_min, y = mean_density)) +
    geom_line(color = sp_col, linewidth = 0.6) +
    labs(title = sprintf("Mean Nuclear Density per Voxel (%d um)", VOXEL_SIZE_UM),
         x = "Time (min)", y = "Mean nuclei / voxel") +
    theme_pipe()

  p_dens_iqr <- ggplot(as.data.frame(dens_stats), aes(x = time_min)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = sp_col, alpha = 0.3) +
    geom_line(aes(y = median_density), color = sp_col, linewidth = 0.7) +
    labs(title = "Voxel Density Over Time",
         subtitle = "Median +/- IQR",
         x = "Time (min)", y = "Nuclei per voxel") +
    theme_pipe()

  p_dens_max <- ggplot(as.data.frame(dens_stats), aes(x = time_min, y = max_density)) +
    geom_line(color = sp_col, alpha = 0.5, linewidth = 0.5) +
    geom_smooth(method = "loess", span = 0.3, se = FALSE,
                color = sp_col, linewidth = 0.8) +
    labs(title = "Max Voxel Density (Hotspot)",
         x = "Time (min)", y = "Max nuclei in single voxel") +
    theme_pipe()

  # IN_ROI density over time
  has_roi <- "IN_ROI" %in% names(spots) && any(spots$IN_ROI == TRUE, na.rm = TRUE)
  if (has_roi) {
    roi_density <- dt[, .(
      n_in_roi     = sum(IN_ROI == TRUE, na.rm = TRUE),
      n_outside    = sum(IN_ROI == FALSE | is.na(IN_ROI)),
      total        = .N
    ), by = FRAME]
    roi_density[, `:=`(frac_roi = n_in_roi / total,
                        time_min = frame_to_min(FRAME))]

    p_roi_frac <- ggplot(as.data.frame(roi_density),
                         aes(x = time_min, y = frac_roi)) +
      geom_line(color = sp_col, linewidth = 0.6) +
      scale_y_continuous(labels = percent) +
      labs(title = "Fraction of Nuclei in ROI Over Time",
           x = "Time (min)", y = "Fraction in ROI") +
      theme_pipe()

    p_roi_count <- ggplot(as.data.frame(roi_density), aes(x = time_min)) +
      geom_line(aes(y = n_in_roi), color = sp_col, linewidth = 0.6) +
      geom_line(aes(y = n_outside), color = "grey50", linewidth = 0.4,
                linetype = "dashed") +
      labs(title = "Nuclei Count: ROI vs. Outside",
           subtitle = "Solid = in ROI, Dashed = outside",
           x = "Time (min)", y = "Number of nuclei") +
      theme_pipe()
  } else {
    p_roi_frac  <- plot_spacer()
    p_roi_count <- plot_spacer()
  }

  p04 <- (p_dens_mean | p_dens_iqr) / (p_dens_max | p_roi_frac) / p_roi_count +
    plot_annotation(
      title = sprintf("%s -- 04 Nuclear Density", sp$name),
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

  save_pdf(p04, sprintf("%s_04_density.pdf", prefix), width = 14, height = 14)


  # =========================================================================
  # S05: VELOCITY ANALYSIS
  # =========================================================================

  cat(sprintf("\n=== %s S05: VELOCITY ANALYSIS ===\n", sp$short))

  # Radial velocity distribution
  p_vr_hist <- spots %>%
    filter(!is.na(RADIAL_VELOCITY_SMOOTH)) %>%
    ggplot(aes(x = RADIAL_VELOCITY_SMOOTH)) +
    geom_histogram(bins = 100, fill = sp_col, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    coord_cartesian(xlim = c(-5, 5)) +
    labs(title = "Radial Velocity Distribution",
         subtitle = "Negative = inward, Positive = outward",
         x = "Smoothed radial velocity (um/frame)", y = "Count") +
    theme_pipe()

  # Instantaneous 3D speed distribution
  p_speed_hist <- spots_v %>%
    filter(!is.na(inst_speed_um_min), inst_speed_um_min < 20) %>%
    ggplot(aes(x = inst_speed_um_min)) +
    geom_histogram(bins = 100, fill = sp_col, alpha = 0.7) +
    labs(title = "Instantaneous 3D Speed Distribution",
         subtitle = "Step-to-step displacement / dt",
         x = "Speed (um/min)", y = "Count") +
    theme_pipe()

  # Tangential vs radial speed scatter
  p_vr_vt <- spots_v %>%
    filter(!is.na(tangential_speed), !is.na(RADIAL_VELOCITY)) %>%
    sample_n(min(100000, n())) %>%
    ggplot(aes(x = abs(RADIAL_VELOCITY), y = tangential_speed)) +
    geom_bin2d(bins = 80) +
    scale_fill_viridis_c(option = "magma", trans = "log10") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "white") +
    coord_cartesian(xlim = c(0, 5), ylim = c(0, 5)) +
    labs(title = "Tangential vs. |Radial| Speed",
         subtitle = "Above diagonal = tangential-dominated movement",
         x = "|Radial velocity| (um/frame)", y = "Tangential speed (um/frame)",
         fill = "Count") +
    theme_pipe()

  # Radial velocity over time
  vel_time <- spots %>%
    filter(!is.na(RADIAL_VELOCITY_SMOOTH)) %>%
    group_by(FRAME) %>%
    summarise(
      mean_vr   = mean(RADIAL_VELOCITY_SMOOTH),
      median_vr = median(RADIAL_VELOCITY_SMOOTH),
      q25_vr    = quantile(RADIAL_VELOCITY_SMOOTH, 0.25),
      q75_vr    = quantile(RADIAL_VELOCITY_SMOOTH, 0.75),
      .groups   = "drop"
    ) %>%
    mutate(time_min = frame_to_min(FRAME))

  p_vr_time <- ggplot(vel_time, aes(x = time_min)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(ymin = q25_vr, ymax = q75_vr), fill = sp_col, alpha = 0.2) +
    geom_line(aes(y = median_vr), color = sp_col, linewidth = 0.7) +
    labs(title = "Radial Velocity Over Time",
         subtitle = "Median +/- IQR",
         x = "Time (min)", y = "Radial velocity (um/frame)") +
    theme_pipe()

  # Instantaneous speed over time (binned)
  speed_time <- spots_v %>%
    filter(!is.na(inst_speed_um_min)) %>%
    mutate(time_bin = floor(frame_to_min(FRAME) / 5) * 5) %>%
    group_by(time_bin) %>%
    summarise(mean_speed = mean(inst_speed_um_min),
              median_speed = median(inst_speed_um_min),
              .groups = "drop")

  p_speed_time <- ggplot(speed_time, aes(x = time_bin, y = median_speed)) +
    geom_line(color = sp_col, linewidth = 0.7) +
    geom_point(color = sp_col, size = 1) +
    labs(title = "Instantaneous Speed Over Time (5-min bins)",
         x = "Time (min)", y = "Median speed (um/min)") +
    theme_pipe()

  # Velocity by depth layer
  p_vr_depth <- add_depth_bin(spots) %>%
    filter(!is.na(RADIAL_VELOCITY_SMOOTH)) %>%
    ggplot(aes(x = depth_bin, y = RADIAL_VELOCITY_SMOOTH)) +
    geom_boxplot(fill = sp_col, alpha = 0.5, outlier.alpha = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    coord_cartesian(ylim = c(-3, 3)) +
    labs(title = "Radial Velocity by Depth Layer",
         x = "Depth layer", y = "Radial velocity (um/frame)") +
    theme_pipe() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # Per-track velocity distributions
  p_track_vr <- tracks %>%
    filter(!is.na(TRACK_MEDIAN_RADIAL_VEL)) %>%
    ggplot(aes(x = TRACK_MEDIAN_RADIAL_VEL)) +
    geom_histogram(bins = 80, fill = sp_col, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    coord_cartesian(xlim = c(-3, 3)) +
    labs(title = "Per-Track Median Radial Velocity",
         x = "Median radial velocity (um/frame)", y = "Count") +
    theme_pipe()

  # Radial velocity x angle heatmap
  p_vr_angle_heat <- spots %>%
    filter(!is.na(RADIAL_VELOCITY_SMOOTH), !is.na(THETA_DEG)) %>%
    ggplot(aes(x = THETA_DEG, y = RADIAL_VELOCITY_SMOOTH)) +
    geom_bin2d(bins = c(50, 60)) +
    scale_fill_viridis_c(option = "inferno", trans = "log10") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "white") +
    coord_cartesian(ylim = c(-4, 4)) +
    labs(title = "Radial Velocity vs. Polar Angle",
         x = "Theta (polar angle, deg)", y = "Radial velocity (um/frame)", fill = "Count") +
    theme_pipe()

  p05 <- (p_vr_hist | p_speed_hist | p_vr_vt) /
         (p_vr_time | p_speed_time) /
         (p_vr_depth | p_track_vr) /
         p_vr_angle_heat +
    plot_annotation(
      title = sprintf("%s -- 05 Velocity Analysis", sp$name),
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

  save_pdf(p05, sprintf("%s_05_velocity.pdf", prefix), width = 18, height = 22)


  # =========================================================================
  # S06: FLOW FIELDS (VELOCITY VECTOR FIELDS)
  # =========================================================================

  cat(sprintf("\n=== %s S06: FLOW FIELDS ===\n", sp$short))

  # --- Theta-Phi velocity field ---
  flow_tp <- spots_v %>%
    filter(!is.na(THETA_DEG), !is.na(PHI_DEG), !is.na(RADIAL_VELOCITY)) %>%
    mutate(theta_bin = floor(THETA_DEG / FLOW_BIN_THETA) * FLOW_BIN_THETA + FLOW_BIN_THETA / 2,
           phi_bin   = floor(PHI_DEG / FLOW_BIN_PHI) * FLOW_BIN_PHI + FLOW_BIN_PHI / 2) %>%
    group_by(theta_bin, phi_bin) %>%
    summarise(
      mean_vr    = mean(RADIAL_VELOCITY_SMOOTH, na.rm = TRUE),
      mean_vt    = mean(tangential_speed, na.rm = TRUE),
      median_vr  = median(RADIAL_VELOCITY_SMOOTH, na.rm = TRUE),
      mean_speed = mean(total_speed_frame, na.rm = TRUE),
      inward_frac = mean(RADIAL_VELOCITY < 0, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(n >= FLOW_MIN_N)

  # Radial velocity heatmap on theta-phi grid
  p_vr_tp_heat <- ggplot(flow_tp, aes(x = phi_bin, y = theta_bin, fill = mean_vr)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1.5, 1.5),
                         oob = squish) +
    labs(title = "Mean Radial Velocity on Theta-Phi Grid",
         subtitle = "Blue = inward, Red = outward",
         x = "Phi (azimuthal, deg)", y = "Theta (polar, deg)",
         fill = "V_r (um/frame)") +
    theme_pipe()

  # Inward fraction heatmap
  p_inward_tp <- ggplot(flow_tp, aes(x = phi_bin, y = theta_bin, fill = inward_frac)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", limits = c(0.3, 0.7), oob = squish,
                         labels = percent) +
    labs(title = "Inward Fraction on Theta-Phi Grid",
         subtitle = ">50% = predominantly inward flow",
         x = "Phi (azimuthal, deg)", y = "Theta (polar, deg)",
         fill = "Inward frac") +
    theme_pipe()

  # Speed magnitude on theta-phi grid
  p_speed_tp <- ggplot(flow_tp, aes(x = phi_bin, y = theta_bin, fill = mean_speed)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma") +
    labs(title = "Mean Speed Magnitude on Theta-Phi Grid",
         x = "Phi (azimuthal, deg)", y = "Theta (polar, deg)",
         fill = "Speed (um/frame)") +
    theme_pipe()

  # --- XY velocity field (arrows) ---
  flow_xy <- spots_v %>%
    filter(!is.na(vx), !is.na(vy)) %>%
    mutate(x_bin = floor(POSITION_X / FLOW_BIN_XY) * FLOW_BIN_XY + FLOW_BIN_XY / 2,
           y_bin = floor(POSITION_Y / FLOW_BIN_XY) * FLOW_BIN_XY + FLOW_BIN_XY / 2) %>%
    group_by(x_bin, y_bin) %>%
    summarise(
      mean_vx    = mean(vx, na.rm = TRUE),
      mean_vy    = mean(vy, na.rm = TRUE),
      mean_speed = mean(total_speed_frame, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(n >= FLOW_MIN_N)

  p_flow_xy <- ggplot(flow_xy) +
    geom_segment(aes(x = x_bin, y = y_bin,
                     xend = x_bin + mean_vx * 10,
                     yend = y_bin + mean_vy * 10,
                     color = mean_speed),
                 arrow = arrow(length = unit(0.08, "cm")),
                 linewidth = 0.4) +
    scale_color_viridis_c(option = "plasma") +
    labs(title = "XY Velocity Vector Field",
         subtitle = "Arrows: mean displacement direction x 10",
         x = "X (um)", y = "Y (um)", color = "Speed (um/frame)") +
    theme_pipe()

  # XY speed heatmap
  p_speed_xy <- flow_xy %>%
    ggplot(aes(x = x_bin, y = y_bin, fill = mean_speed)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma") +
    labs(title = "XY Speed Heatmap",
         x = "X (um)", y = "Y (um)", fill = "Speed") +
    theme_pipe()

  # Temporal evolution of flow: radial velocity at time windows
  time_range <- range(spots$FRAME)
  n_windows <- 4
  window_size <- ceiling(diff(time_range) / n_windows)
  flow_windows <- spots_v %>%
    filter(!is.na(THETA_DEG), !is.na(RADIAL_VELOCITY_SMOOTH)) %>%
    mutate(
      win_idx = pmin(floor((FRAME - time_range[1]) / window_size), n_windows - 1),
      time_window = paste0(
        sprintf("%.0f", frame_to_min(win_idx * window_size + time_range[1])),
        "-",
        sprintf("%.0f min", frame_to_min((win_idx + 1) * window_size + time_range[1]))
      )
    ) %>%
    mutate(theta_bin = floor(THETA_DEG / FLOW_BIN_THETA) * FLOW_BIN_THETA + FLOW_BIN_THETA / 2,
           phi_bin   = floor(PHI_DEG / FLOW_BIN_PHI) * FLOW_BIN_PHI + FLOW_BIN_PHI / 2) %>%
    group_by(theta_bin, phi_bin, time_window) %>%
    summarise(mean_vr = mean(RADIAL_VELOCITY_SMOOTH, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    filter(n >= 10)

  p_flow_time <- ggplot(flow_windows, aes(x = phi_bin, y = theta_bin, fill = mean_vr)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1.5, 1.5), oob = squish) +
    facet_wrap(~time_window, ncol = 2) +
    labs(title = "Radial Velocity Field -- Temporal Evolution",
         x = "Phi (deg)", y = "Theta (deg)", fill = "V_r") +
    theme_pipe()

  p06 <- (p_vr_tp_heat | p_inward_tp) / (p_speed_tp | p_flow_xy) / (p_speed_xy | plot_spacer()) / p_flow_time +
    plot_annotation(
      title = sprintf("%s -- 06 Flow Fields", sp$name),
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

  save_pdf(p06, sprintf("%s_06_flow_fields.pdf", prefix), width = 16, height = 24)


  # =========================================================================
  # S07: PHYSICAL PROPERTIES (MSD, CONFINEMENT, DISPLACEMENT)
  # =========================================================================

  cat(sprintf("\n=== %s S07: PHYSICAL PROPERTIES ===\n", sp$short))

  # --- MSD analysis ---
  cat("    Computing MSD ...\n")
  msd_df <- compute_msd_ensemble(spots, max_lag = MSD_MAX_LAG, n_sample = MSD_N_TRACKS)

  # Fit power law: MSD = D * t^alpha on log-log
  msd_fit_data <- msd_df %>% filter(lag_min <= 10, mean_msd > 0)
  if (nrow(msd_fit_data) >= 3) {
    msd_fit <- lm(log(mean_msd) ~ log(lag_min), data = msd_fit_data)
    alpha <- coef(msd_fit)[2]
    D_apparent <- exp(coef(msd_fit)[1]) / 6  # 3D diffusion: MSD = 6D*t^alpha
    alpha_label <- ifelse(alpha < 0.9, "Subdiffusive (confined)",
                   ifelse(alpha < 1.1, "Normal diffusion",
                          "Superdiffusive (directed)"))
  } else {
    alpha <- NA; D_apparent <- NA; alpha_label <- "N/A"
  }
  cat(sprintf("    MSD alpha = %.3f -> %s, D_app = %.3f um^2/min\n",
              alpha, alpha_label, D_apparent))

  p_msd <- ggplot(msd_df, aes(x = lag_min, y = mean_msd)) +
    geom_point(color = sp_col, size = 1.5) +
    geom_line(color = sp_col, linewidth = 0.5) +
    { if (!is.na(alpha)) geom_smooth(method = "lm", formula = y ~ x,
                                      se = FALSE, color = "black",
                                      linetype = "dashed", linewidth = 0.6) } +
    scale_x_log10() + scale_y_log10() +
    labs(title = "Mean Squared Displacement (MSD)",
         subtitle = sprintf("alpha = %.3f (%s), D_app = %.3f um^2/min",
                            alpha, alpha_label, D_apparent),
         x = "Lag time (min, log)", y = "MSD (um^2, log)") +
    theme_pipe()

  p_msd_linear <- ggplot(msd_df, aes(x = lag_min, y = mean_msd)) +
    geom_point(color = sp_col, size = 1.5) +
    geom_line(color = sp_col, linewidth = 0.5) +
    labs(title = "MSD (Linear Scale)",
         subtitle = "Curvature indicates anomalous diffusion",
         x = "Lag time (min)", y = "MSD (um^2)") +
    theme_pipe()

  # --- Net radial displacement ---
  p_netdisp <- tracks %>%
    filter(!is.na(TRACK_NET_RADIAL_DISP)) %>%
    ggplot(aes(x = TRACK_NET_RADIAL_DISP)) +
    geom_histogram(bins = 80, fill = sp_col, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    coord_cartesian(xlim = quantile(tracks$TRACK_NET_RADIAL_DISP,
                                    c(0.01, 0.99), na.rm = TRUE)) +
    labs(title = "Net Radial Displacement per Track",
         subtitle = "Negative = net inward, Positive = net outward",
         x = "Net radial displacement (um)", y = "Count") +
    theme_pipe()

  # --- Displacement vs track length colored by depth ---
  p_disp_len <- tracks %>%
    filter(!is.na(TRACK_NET_RADIAL_DISP)) %>%
    ggplot(aes(x = n_spots, y = TRACK_NET_RADIAL_DISP, color = mean_SPHERICAL_DEPTH)) +
    geom_point(alpha = 0.08, size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    scale_color_viridis_c(option = "mako", direction = -1) +
    scale_x_log10() +
    labs(title = "Net Displacement vs. Track Length",
         subtitle = "Coloured by mean spherical depth",
         x = "Track length (spots, log)", y = "Net radial displacement (um)",
         color = "Depth") +
    theme_pipe()

  # --- Cumulative displacement distribution ---
  cum_disp <- tracks %>%
    filter(!is.na(TRACK_NET_RADIAL_DISP)) %>%
    arrange(TRACK_NET_RADIAL_DISP) %>%
    mutate(frac_below = row_number() / n())

  p_cum_disp <- ggplot(cum_disp, aes(x = TRACK_NET_RADIAL_DISP, y = frac_below)) +
    geom_line(color = sp_col, linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    scale_y_continuous(labels = percent) +
    labs(title = "Cumulative Distribution of Net Displacement",
         subtitle = "Left half = net inward tracks",
         x = "Net radial displacement (um)", y = "Cumulative fraction") +
    theme_pipe()

  # --- Instantaneous speed vs depth scatter ---
  p_speed_depth <- spots_v %>%
    filter(!is.na(inst_speed_um_min), !is.na(SPHERICAL_DEPTH),
           inst_speed_um_min < 15) %>%
    sample_n(min(100000, n())) %>%
    ggplot(aes(x = SPHERICAL_DEPTH, y = inst_speed_um_min)) +
    geom_bin2d(bins = c(50, 50)) +
    scale_fill_viridis_c(option = "magma", trans = "log10") +
    labs(title = "Speed vs. Depth",
         x = "Spherical depth (um)", y = "Inst. speed (um/min)",
         fill = "Count") +
    theme_pipe()

  p07 <- (p_msd | p_msd_linear) / (p_netdisp | p_cum_disp) / (p_disp_len | p_speed_depth) +
    plot_annotation(
      title = sprintf("%s -- 07 Physical Properties", sp$name),
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

  save_pdf(p07, sprintf("%s_07_physical_properties.pdf", prefix), width = 16, height = 16)


  # =========================================================================
  # S08: DIRECTIONALITY
  # =========================================================================

  cat(sprintf("\n=== %s S08: DIRECTIONALITY ===\n", sp$short))

  # Inward fraction distribution per track
  p_inward_dist <- tracks %>%
    filter(!is.na(TRACK_INWARD_FRACTION)) %>%
    ggplot(aes(x = TRACK_INWARD_FRACTION)) +
    geom_histogram(bins = 50, fill = sp_col, alpha = 0.7) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey30") +
    labs(title = "Inward Fraction per Track",
         subtitle = ">0.5 = predominantly inward",
         x = "Inward fraction", y = "Count") +
    theme_pipe()

  # Sustained inward movement
  p_sustained <- tracks %>%
    filter(!is.na(TRACK_SUSTAINED_INWARD), TRACK_SUSTAINED_INWARD > 0) %>%
    ggplot(aes(x = TRACK_SUSTAINED_INWARD)) +
    geom_histogram(bins = 50, fill = sp_col, alpha = 0.7) +
    labs(title = "Max Consecutive Inward Frames",
         x = "Consecutive inward frames", y = "Count") +
    theme_pipe()

  # Turning angle distribution
  p_turning <- spots_v %>%
    filter(!is.na(turning_angle)) %>%
    ggplot(aes(x = turning_angle)) +
    geom_histogram(bins = 60, fill = sp_col, alpha = 0.7) +
    geom_vline(xintercept = 90, linetype = "dashed", color = "grey30") +
    labs(title = "Turning Angle Distribution",
         subtitle = "90 deg = random, <90 = persistent, >90 = reversing",
         x = "Turning angle (deg)", y = "Count") +
    theme_pipe()

  # Turn angle vs speed
  p_turn_speed <- spots_v %>%
    filter(!is.na(turning_angle), !is.na(inst_speed_um_min),
           inst_speed_um_min < 15) %>%
    sample_n(min(100000, n())) %>%
    ggplot(aes(x = turning_angle, y = inst_speed_um_min)) +
    geom_bin2d(bins = c(60, 50)) +
    scale_fill_viridis_c(option = "magma", trans = "log10") +
    labs(title = "Turning Angle vs. Instantaneous Speed",
         x = "Turning angle (deg)", y = "Speed (um/min)", fill = "Count") +
    theme_pipe()

  # Rose diagram of XY direction
  p_rose <- spots_v %>%
    filter(!is.na(direction_xy)) %>%
    ggplot(aes(x = direction_xy)) +
    geom_histogram(bins = 36, fill = sp_col, alpha = 0.7) +
    coord_polar(start = 0) +
    scale_x_continuous(breaks = seq(-180, 135, 45)) +
    labs(title = "XY Direction Rose Diagram",
         subtitle = "Uniform = isotropic, peaks = preferred axes",
         x = "Direction (deg)", y = "Count") +
    theme_pipe()

  # Mean turning angle by depth
  turn_depth <- spots_v %>%
    filter(!is.na(turning_angle)) %>%
    add_depth_bin() %>%
    group_by(depth_bin) %>%
    summarise(mean_angle = mean(turning_angle),
              median_angle = median(turning_angle),
              n = n(), .groups = "drop")

  p_turn_depth <- ggplot(turn_depth, aes(x = depth_bin, y = median_angle)) +
    geom_col(fill = sp_col, alpha = 0.7) +
    geom_hline(yintercept = 90, linetype = "dashed", color = "grey30") +
    labs(title = "Median Turning Angle by Depth",
         x = "Depth layer", y = "Median angle (deg)") +
    theme_pipe() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  p08 <- (p_inward_dist | p_sustained | p_turning) /
         (p_turn_speed  | p_rose      | p_turn_depth) +
    plot_annotation(
      title = sprintf("%s -- 08 Directionality", sp$name),
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

  save_pdf(p08, sprintf("%s_08_directionality.pdf", prefix), width = 18, height = 12)


  # =========================================================================
  # S09: INTENSITY & SIGNAL QUALITY
  # =========================================================================

  cat(sprintf("\n=== %s S09: INTENSITY QC ===\n", sp$short))

  # Intensity over time
  int_time <- spots %>%
    filter(!is.na(MEAN_INTENSITY_CH1)) %>%
    group_by(FRAME) %>%
    summarise(median_int = median(MEAN_INTENSITY_CH1),
              mean_snr = mean(SNR_CH1, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(time_min = frame_to_min(FRAME))

  p_int <- ggplot(int_time, aes(x = time_min, y = median_int)) +
    geom_line(color = sp_col, linewidth = 0.6) +
    labs(title = "Median Spot Intensity Over Time",
         x = "Time (min)", y = "Median intensity (CH1)") +
    theme_pipe()

  p_snr <- ggplot(int_time, aes(x = time_min, y = mean_snr)) +
    geom_line(color = sp_col, linewidth = 0.6) +
    labs(title = "Mean SNR Over Time",
         x = "Time (min)", y = "Mean SNR (CH1)") +
    theme_pipe()

  # Intensity vs depth
  p_int_depth <- spots %>%
    filter(!is.na(MEAN_INTENSITY_CH1)) %>%
    ggplot(aes(x = SPHERICAL_DEPTH, y = MEAN_INTENSITY_CH1)) +
    geom_bin2d(bins = c(50, 50)) +
    scale_fill_viridis_c(option = "magma", trans = "log10") +
    labs(title = "Intensity vs. Depth",
         subtitle = "Deeper nuclei likely have lower signal (scattering)",
         x = "Spherical depth (um)", y = "Mean intensity (CH1)", fill = "Count") +
    theme_pipe()

  # Quality score
  p_quality <- spots %>%
    ggplot(aes(x = QUALITY)) +
    geom_histogram(bins = 80, fill = sp_col, alpha = 0.7) +
    scale_x_log10(labels = comma) +
    labs(title = "Spot Quality Score Distribution",
         x = "Quality (log)", y = "Count") +
    theme_pipe()

  p09 <- (p_int | p_snr) / (p_int_depth | p_quality) +
    plot_annotation(
      title = sprintf("%s -- 09 Intensity & Signal QC", sp$name),
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

  save_pdf(p09, sprintf("%s_09_intensity_qc.pdf", prefix), width = 14, height = 10)


  # =========================================================================
  # S10: COMPREHENSIVE INGRESSION ANALYSIS (validated species only)
  # =========================================================================

  if (sp$has_ingression) {

    cat(sprintf("\n=== %s S10: INGRESSION ANALYSIS ===\n", sp$short))

    ingr_tracks <- tracks %>% filter(INGRESSING == TRUE)
    non_ingr_tracks <- tracks %>% filter(INGRESSING == FALSE | is.na(INGRESSING))
    n_ingr <- nrow(ingr_tracks)
    n_total <- nrow(tracks)

    cat(sprintf("    Ingressing: %d / %d tracks (%.1f%%)\n",
                n_ingr, n_total, 100 * n_ingr / n_total))

    ingr_ids <- ingr_tracks$TRACK_ID
    spots_ingr_flag <- spots %>%
      mutate(is_ingressing = TRACK_ID %in% ingr_ids)

    onset_frames <- ingr_tracks %>%
      filter(!is.na(INGRESSION_ONSET_FRAME)) %>%
      select(TRACK_ID, onset = INGRESSION_ONSET_FRAME)

    # --- 10a: Onset timing ---
    p_onset <- ingr_tracks %>%
      filter(!is.na(INGRESSION_ONSET_FRAME)) %>%
      mutate(onset_min = frame_to_min(INGRESSION_ONSET_FRAME)) %>%
      ggplot(aes(x = onset_min)) +
      geom_histogram(bins = 40, fill = "#E41A1C", alpha = 0.7) +
      labs(title = "Ingression Onset Timing",
           subtitle = sprintf("%d ingressing tracks, %d with known onset",
                              n_ingr, nrow(onset_frames)),
           x = "Onset time (min)", y = "Count") +
      theme_pipe()

    # --- 10b: Depth at onset ---
    onset_spots <- spots_ingr_flag %>%
      filter(is_ingressing, !is.na(INGRESSION_ONSET_FRAME)) %>%
      group_by(TRACK_ID) %>%
      arrange(FRAME) %>%
      filter(FRAME == INGRESSION_ONSET_FRAME[1]) %>%
      slice(1) %>%
      ungroup()

    p_onset_depth <- onset_spots %>%
      ggplot(aes(x = SPHERICAL_DEPTH)) +
      geom_histogram(bins = 30, fill = "#E41A1C", alpha = 0.7) +
      labs(title = "Spherical Depth at Ingression Onset",
           x = "Spherical depth (um)", y = "Count") +
      theme_pipe()

    # --- 10c: Angular position of ingression onset ---
    p_onset_angular <- onset_spots %>%
      filter(!is.na(THETA_DEG), !is.na(PHI_DEG)) %>%
      ggplot(aes(x = PHI_DEG, y = THETA_DEG)) +
      geom_bin2d(bins = c(30, 20)) +
      scale_fill_viridis_c(option = "inferno") +
      labs(title = "Angular Position at Ingression Onset",
           subtitle = "Where on the embryo surface do cells ingress?",
           x = "Phi (deg)", y = "Theta (deg)", fill = "Count") +
      theme_pipe()

    # --- 10d: Spatial distribution ingressing vs non-ingressing (theta-depth) ---
    p_spatial_ingr <- spots_ingr_flag %>%
      ggplot(aes(x = THETA_DEG, y = SPHERICAL_DEPTH, color = is_ingressing)) +
      geom_point(data = . %>% filter(!is_ingressing),
                 alpha = 0.003, size = 0.1, color = "grey70") +
      geom_point(data = . %>% filter(is_ingressing),
                 alpha = 0.02, size = 0.3, color = "#E41A1C") +
      labs(title = "Ingressing (red) vs. Non-Ingressing (grey)",
           subtitle = "Theta-depth projection",
           x = "Theta (polar, deg)", y = "Spherical depth (um)") +
      theme_pipe() + theme(legend.position = "none")

    # --- 10e: Ingression score distribution ---
    p_score <- tracks %>%
      filter(!is.na(INGRESSION_SCORE)) %>%
      mutate(label = ifelse(INGRESSING == TRUE, "Ingressing", "Non-ingressing")) %>%
      ggplot(aes(x = INGRESSION_SCORE, fill = label)) +
      geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
      scale_fill_manual(values = c("Ingressing" = "#E41A1C",
                                   "Non-ingressing" = "grey60")) +
      labs(title = "Ingression Score Distribution",
           x = "Score (0 = none, 1 = strong)", y = "Count", fill = "") +
      theme_pipe()

    # --- 10f: Ingression rate over time (5-min bins) ---
    ingr_rate <- ingr_tracks %>%
      filter(!is.na(INGRESSION_ONSET_FRAME)) %>%
      mutate(onset_min = frame_to_min(INGRESSION_ONSET_FRAME),
             time_bin = floor(onset_min / 5) * 5) %>%
      group_by(time_bin) %>%
      summarise(n_new = n(), .groups = "drop")

    p_rate_ingr <- ggplot(ingr_rate, aes(x = time_bin, y = n_new)) +
      geom_col(fill = "#E41A1C", alpha = 0.7) +
      labs(title = "Ingression Rate Over Time",
           subtitle = "New onsets per 5-min window",
           x = "Time (min)", y = "New ingressing tracks") +
      theme_pipe()

    # Cumulative ingression over time
    ingr_cum <- ingr_tracks %>%
      filter(!is.na(INGRESSION_ONSET_FRAME)) %>%
      arrange(INGRESSION_ONSET_FRAME) %>%
      mutate(cum_n = row_number(),
             onset_min = frame_to_min(INGRESSION_ONSET_FRAME))

    p_cum_ingr <- ggplot(ingr_cum, aes(x = onset_min, y = cum_n)) +
      geom_line(color = "#E41A1C", linewidth = 0.8) +
      labs(title = "Cumulative Ingression Onsets",
           x = "Time (min)", y = "Cumulative ingressing tracks") +
      theme_pipe()

    # --- 10g: Velocity pseudo-trajectory aligned to onset ---
    aligned_vel <- spots %>%
      filter(TRACK_ID %in% onset_frames$TRACK_ID,
             !is.na(RADIAL_VELOCITY_SMOOTH)) %>%
      inner_join(onset_frames, by = "TRACK_ID") %>%
      mutate(rel_frame = FRAME - onset)

    vel_profile <- aligned_vel %>%
      group_by(rel_frame) %>%
      summarise(
        mean_vr   = mean(RADIAL_VELOCITY_SMOOTH),
        median_vr = median(RADIAL_VELOCITY_SMOOTH),
        se_vr     = sd(RADIAL_VELOCITY_SMOOTH) / sqrt(n()),
        q25_vr    = quantile(RADIAL_VELOCITY_SMOOTH, 0.25),
        q75_vr    = quantile(RADIAL_VELOCITY_SMOOTH, 0.75),
        n         = n(),
        .groups   = "drop"
      ) %>%
      filter(n >= 5, abs(rel_frame) <= 30)

    p_vel_pseudo <- ggplot(vel_profile,
                           aes(x = rel_frame * FRAME_INTERVAL_MIN, y = mean_vr)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#E41A1C",
                 linewidth = 0.8) +
      geom_ribbon(aes(ymin = q25_vr, ymax = q75_vr),
                  fill = "#E41A1C", alpha = 0.15) +
      geom_ribbon(aes(ymin = mean_vr - se_vr, ymax = mean_vr + se_vr),
                  fill = "#E41A1C", alpha = 0.3) +
      geom_line(linewidth = 0.8, color = "#E41A1C") +
      annotate("text", x = 0.3, y = max(vel_profile$mean_vr, na.rm = TRUE) * 0.8,
               label = "Onset", color = "#E41A1C", hjust = 0, fontface = "italic") +
      labs(title = "Radial Velocity -- Pseudo-Trajectory",
           subtitle = "Aligned to ingression onset (mean +/- SE, ribbon = IQR)",
           x = "Time relative to onset (min)", y = "Radial velocity (um/frame)") +
      theme_pipe()

    # --- 10h: Depth pseudo-trajectory aligned to onset ---
    depth_traj <- spots %>%
      filter(TRACK_ID %in% onset_frames$TRACK_ID) %>%
      inner_join(onset_frames, by = "TRACK_ID") %>%
      mutate(rel_frame = FRAME - onset)

    depth_profile <- depth_traj %>%
      group_by(rel_frame) %>%
      summarise(
        mean_depth   = mean(SPHERICAL_DEPTH, na.rm = TRUE),
        median_depth = median(SPHERICAL_DEPTH, na.rm = TRUE),
        q25_depth    = quantile(SPHERICAL_DEPTH, 0.25, na.rm = TRUE),
        q75_depth    = quantile(SPHERICAL_DEPTH, 0.75, na.rm = TRUE),
        n = n(), .groups = "drop"
      ) %>%
      filter(n >= 5, abs(rel_frame) <= 30)

    p_depth_pseudo <- ggplot(depth_profile,
                             aes(x = rel_frame * FRAME_INTERVAL_MIN)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#E41A1C",
                 linewidth = 0.8) +
      geom_ribbon(aes(ymin = q25_depth, ymax = q75_depth),
                  fill = "#E41A1C", alpha = 0.15) +
      geom_line(aes(y = median_depth), color = "#E41A1C", linewidth = 0.8) +
      labs(title = "Depth -- Pseudo-Trajectory",
           subtitle = "Median +/- IQR, aligned to onset",
           x = "Time relative to onset (min)", y = "Spherical depth (um)") +
      theme_pipe()

    # --- 10i: Radial distance pseudo-trajectory ---
    radial_profile <- depth_traj %>%
      group_by(rel_frame) %>%
      summarise(
        mean_rdist   = mean(RADIAL_DIST_TO_CENTER, na.rm = TRUE),
        median_rdist = median(RADIAL_DIST_TO_CENTER, na.rm = TRUE),
        q25_rdist    = quantile(RADIAL_DIST_TO_CENTER, 0.25, na.rm = TRUE),
        q75_rdist    = quantile(RADIAL_DIST_TO_CENTER, 0.75, na.rm = TRUE),
        n = n(), .groups = "drop"
      ) %>%
      filter(n >= 5, abs(rel_frame) <= 30)

    p_rdist_pseudo <- ggplot(radial_profile,
                             aes(x = rel_frame * FRAME_INTERVAL_MIN)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#E41A1C",
                 linewidth = 0.8) +
      geom_ribbon(aes(ymin = q25_rdist, ymax = q75_rdist),
                  fill = "#E41A1C", alpha = 0.15) +
      geom_line(aes(y = median_rdist), color = "#E41A1C", linewidth = 0.8) +
      geom_hline(yintercept = R_sphere, linetype = "dotted", color = "grey40") +
      annotate("text", x = min(radial_profile$rel_frame) * FRAME_INTERVAL_MIN,
               y = R_sphere * 1.01, label = "R_sphere", hjust = 0,
               color = "grey40", size = 3, fontface = "italic") +
      labs(title = "Radial Distance -- Pseudo-Trajectory",
           subtitle = "Distance from sphere centre; dashed = sphere surface",
           x = "Time relative to onset (min)",
           y = "Radial distance to centre (um)") +
      theme_pipe()

    # --- 10j: Speed pseudo-trajectory (absolute inst. speed around onset) ---
    aligned_speed <- d$spots_vel %>%
      filter(TRACK_ID %in% onset_frames$TRACK_ID,
             !is.na(inst_speed_um_min)) %>%
      inner_join(onset_frames, by = "TRACK_ID") %>%
      mutate(rel_frame = FRAME - onset)

    speed_profile <- aligned_speed %>%
      group_by(rel_frame) %>%
      summarise(
        mean_speed   = mean(inst_speed_um_min),
        median_speed = median(inst_speed_um_min),
        q25          = quantile(inst_speed_um_min, 0.25),
        q75          = quantile(inst_speed_um_min, 0.75),
        n = n(), .groups = "drop"
      ) %>%
      filter(n >= 5, abs(rel_frame) <= 30)

    p_speed_pseudo <- ggplot(speed_profile,
                             aes(x = rel_frame * FRAME_INTERVAL_MIN)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#E41A1C",
                 linewidth = 0.8) +
      geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#E41A1C", alpha = 0.15) +
      geom_line(aes(y = median_speed), color = "#E41A1C", linewidth = 0.8) +
      labs(title = "3D Speed -- Pseudo-Trajectory",
           subtitle = "Median +/- IQR, aligned to onset",
           x = "Time relative to onset (min)", y = "Speed (um/min)") +
      theme_pipe()

    # --- 10k: Pre-ingression signatures ---
    # Compare velocity 5 frames before onset to all non-ingressing at same time
    pre_onset_window <- 5  # frames before onset
    pre_ingr <- aligned_vel %>%
      filter(rel_frame >= -pre_onset_window, rel_frame < 0) %>%
      group_by(TRACK_ID) %>%
      summarise(pre_mean_vr = mean(RADIAL_VELOCITY_SMOOTH),
                .groups = "drop") %>%
      mutate(group = "Pre-ingression")

    # Non-ingressing tracks: average velocity
    non_ingr_vel <- spots %>%
      filter(!TRACK_ID %in% ingr_ids, !is.na(RADIAL_VELOCITY_SMOOTH)) %>%
      group_by(TRACK_ID) %>%
      summarise(pre_mean_vr = mean(RADIAL_VELOCITY_SMOOTH),
                .groups = "drop") %>%
      mutate(group = "Non-ingressing")

    pre_comp <- bind_rows(pre_ingr, non_ingr_vel)

    p_pre_sig <- ggplot(pre_comp, aes(x = pre_mean_vr, fill = group)) +
      geom_density(alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_fill_manual(values = c("Pre-ingression" = "#E41A1C",
                                   "Non-ingressing" = "grey60")) +
      coord_cartesian(xlim = c(-3, 3)) +
      labs(title = "Pre-Ingression Velocity Signature",
           subtitle = sprintf("Mean V_r in %d frames before onset vs. all non-ingressing",
                              pre_onset_window),
           x = "Mean radial velocity (um/frame)", y = "Density", fill = "") +
      theme_pipe()

    # --- 10l: Comparison violin plots (ingr vs non-ingr) ---
    comp_data <- tracks %>%
      mutate(group = ifelse(INGRESSING == TRUE, "Ingressing", "Non-ingressing"))

    p_comp_len <- comp_data %>%
      ggplot(aes(x = group, y = n_spots, fill = group)) +
      geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
      scale_fill_manual(values = c("Ingressing" = "#E41A1C",
                                   "Non-ingressing" = "grey60")) +
      labs(title = "Track Length", y = "Spots", x = "") +
      theme_pipe() + theme(legend.position = "none")

    p_comp_depth <- comp_data %>%
      ggplot(aes(x = group, y = mean_SPHERICAL_DEPTH, fill = group)) +
      geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
      scale_fill_manual(values = c("Ingressing" = "#E41A1C",
                                   "Non-ingressing" = "grey60")) +
      labs(title = "Mean Depth", y = "Depth (um)", x = "") +
      theme_pipe() + theme(legend.position = "none")

    p_comp_vel <- comp_data %>%
      filter(!is.na(TRACK_MEDIAN_RADIAL_VEL)) %>%
      ggplot(aes(x = group, y = TRACK_MEDIAN_RADIAL_VEL, fill = group)) +
      geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
      scale_fill_manual(values = c("Ingressing" = "#E41A1C",
                                   "Non-ingressing" = "grey60")) +
      coord_cartesian(ylim = c(-3, 3)) +
      labs(title = "Median V_r", y = "um/frame", x = "") +
      theme_pipe() + theme(legend.position = "none")

    p_comp_inward <- comp_data %>%
      filter(!is.na(TRACK_INWARD_FRACTION)) %>%
      ggplot(aes(x = group, y = TRACK_INWARD_FRACTION, fill = group)) +
      geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
      scale_fill_manual(values = c("Ingressing" = "#E41A1C",
                                   "Non-ingressing" = "grey60")) +
      labs(title = "Inward Fraction", y = "Fraction", x = "") +
      theme_pipe() + theme(legend.position = "none")

    p_comp_theta <- comp_data %>%
      filter(!is.na(mean_THETA_DEG)) %>%
      ggplot(aes(x = group, y = mean_THETA_DEG, fill = group)) +
      geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
      scale_fill_manual(values = c("Ingressing" = "#E41A1C",
                                   "Non-ingressing" = "grey60")) +
      labs(title = "Mean Theta (Polar)", y = "Theta (deg)", x = "") +
      theme_pipe() + theme(legend.position = "none")

    p_comp_score <- comp_data %>%
      filter(!is.na(INGRESSION_SCORE)) %>%
      ggplot(aes(x = group, y = INGRESSION_SCORE, fill = group)) +
      geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
      scale_fill_manual(values = c("Ingressing" = "#E41A1C",
                                   "Non-ingressing" = "grey60")) +
      labs(title = "Ingression Score", y = "Score", x = "") +
      theme_pipe() + theme(legend.position = "none")

    # --- Assemble Section 10 ---
    p10_row1 <- (p_onset | p_onset_depth | p_onset_angular)
    p10_row2 <- (p_spatial_ingr | p_score | p_rate_ingr)
    p10_row3 <- (p_cum_ingr | p_pre_sig)
    p10_row4 <- (p_vel_pseudo | p_depth_pseudo | p_rdist_pseudo)
    p10_row5 <- p_speed_pseudo
    p10_row6 <- (p_comp_len | p_comp_depth | p_comp_vel)
    p10_row7 <- (p_comp_inward | p_comp_theta | p_comp_score)

    p10 <- p10_row1 / p10_row2 / p10_row3 / p10_row4 / p10_row5 / p10_row6 / p10_row7 +
      plot_annotation(
        title = sprintf("%s -- 10 Ingression Analysis", sp$name),
        subtitle = sprintf("%d ingressing / %s total (%.1f%%)",
                           n_ingr, format(n_total, big.mark = ","),
                           100 * n_ingr / n_total),
        theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5, color = "grey50")))

    save_pdf(p10, sprintf("%s_10_ingression.pdf", prefix), width = 18, height = 36)


    # =====================================================================
    # S11: MARGIN / ROI DENSITY ANALYSIS (ingression zone)
    # =====================================================================

    cat(sprintf("\n=== %s S11: MARGIN DENSITY ===\n", sp$short))

    # Density of nuclei in ROI over time
    margin_density <- spots %>%
      mutate(in_roi = as.logical(IN_ROI),
             time_min = frame_to_min(FRAME)) %>%
      group_by(FRAME, time_min, in_roi) %>%
      summarise(n = n(), .groups = "drop") %>%
      pivot_wider(names_from = in_roi, values_from = n,
                  names_prefix = "roi_", values_fill = 0)

    # Handle possible missing columns
    if (!"roi_TRUE" %in% names(margin_density)) margin_density$roi_TRUE <- 0
    if (!"roi_FALSE" %in% names(margin_density)) margin_density$roi_FALSE <- 0

    margin_density <- margin_density %>%
      mutate(total = roi_TRUE + roi_FALSE,
             frac_roi = roi_TRUE / total)

    p_margin_count <- ggplot(margin_density, aes(x = time_min)) +
      geom_line(aes(y = roi_TRUE), color = "#E41A1C", linewidth = 0.7) +
      geom_line(aes(y = roi_FALSE), color = "grey50", linewidth = 0.4,
                linetype = "dashed") +
      labs(title = "Nuclei Count: In ROI (red) vs. Outside (grey)",
           x = "Time (min)", y = "Number of nuclei") +
      theme_pipe()

    p_margin_frac <- ggplot(margin_density, aes(x = time_min, y = frac_roi)) +
      geom_line(color = "#E41A1C", linewidth = 0.7) +
      scale_y_continuous(labels = percent) +
      labs(title = "Fraction of Nuclei in ROI Over Time",
           x = "Time (min)", y = "Fraction in ROI") +
      theme_pipe()

    # Density of ingressing cells in ROI specifically
    ingr_in_roi <- spots_ingr_flag %>%
      mutate(in_roi = as.logical(IN_ROI),
             time_min = frame_to_min(FRAME)) %>%
      filter(in_roi) %>%
      group_by(FRAME, time_min) %>%
      summarise(n_total_roi = n(),
                n_ingr_roi  = sum(is_ingressing),
                .groups = "drop") %>%
      mutate(frac_ingr = n_ingr_roi / n_total_roi)

    p_ingr_roi_frac <- ggplot(ingr_in_roi, aes(x = time_min, y = frac_ingr)) +
      geom_line(color = "#E41A1C", linewidth = 0.5, alpha = 0.5) +
      geom_smooth(method = "loess", span = 0.3, se = TRUE,
                  color = "#E41A1C", fill = "#E41A1C", alpha = 0.2) +
      scale_y_continuous(labels = percent) +
      labs(title = "Fraction of ROI Nuclei That Are Ingressing",
           subtitle = "In the margin zone over time",
           x = "Time (min)", y = "Ingressing / total in ROI") +
      theme_pipe()

    # Depth distribution within ROI: ingressing vs non-ingressing
    p_roi_depth_comp <- spots_ingr_flag %>%
      filter(as.logical(IN_ROI)) %>%
      mutate(label = ifelse(is_ingressing, "Ingressing", "Non-ingressing")) %>%
      ggplot(aes(x = SPHERICAL_DEPTH, fill = label)) +
      geom_density(alpha = 0.5) +
      scale_fill_manual(values = c("Ingressing" = "#E41A1C",
                                   "Non-ingressing" = "grey60")) +
      labs(title = "Depth Distribution Within ROI",
           subtitle = "Ingressing vs. Non-ingressing nuclei in the margin",
           x = "Spherical depth (um)", y = "Density", fill = "") +
      theme_pipe()

    # Theta distribution of all ROI nuclei vs ingressing
    p_roi_theta <- spots_ingr_flag %>%
      filter(as.logical(IN_ROI)) %>%
      mutate(label = ifelse(is_ingressing, "Ingressing", "Non-ingressing")) %>%
      ggplot(aes(x = THETA_DEG, fill = label)) +
      geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
      scale_fill_manual(values = c("Ingressing" = "#E41A1C",
                                   "Non-ingressing" = "grey60")) +
      labs(title = "Theta Distribution Within ROI",
           x = "Theta (polar, deg)", y = "Count", fill = "") +
      theme_pipe()

    # Radial velocity in ROI vs outside
    roi_vel_comp <- spots %>%
      filter(!is.na(RADIAL_VELOCITY_SMOOTH)) %>%
      mutate(location = ifelse(as.logical(IN_ROI), "In ROI", "Outside ROI"))

    p_roi_vel <- ggplot(roi_vel_comp,
                        aes(x = RADIAL_VELOCITY_SMOOTH, fill = location)) +
      geom_density(alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_fill_manual(values = c("In ROI" = "#E41A1C", "Outside ROI" = "grey60")) +
      coord_cartesian(xlim = c(-4, 4)) +
      labs(title = "Radial Velocity: In ROI vs. Outside",
           x = "Radial velocity (um/frame)", y = "Density", fill = "") +
      theme_pipe()

    # Density heatmap on theta-phi: ingressing cell locations
    p_ingr_heat <- spots_ingr_flag %>%
      filter(is_ingressing, !is.na(THETA_DEG), !is.na(PHI_DEG)) %>%
      ggplot(aes(x = PHI_DEG, y = THETA_DEG)) +
      geom_bin2d(bins = c(40, 30)) +
      scale_fill_viridis_c(option = "inferno") +
      labs(title = "Ingressing Cell Density on Theta-Phi Grid",
           x = "Phi (deg)", y = "Theta (deg)", fill = "Count") +
      theme_pipe()

    # Ingressing cell depth over absolute time (not aligned)
    ingr_depth_time <- spots_ingr_flag %>%
      filter(is_ingressing) %>%
      mutate(time_min = frame_to_min(FRAME)) %>%
      group_by(FRAME, time_min) %>%
      summarise(mean_depth = mean(SPHERICAL_DEPTH, na.rm = TRUE),
                median_depth = median(SPHERICAL_DEPTH, na.rm = TRUE),
                q25 = quantile(SPHERICAL_DEPTH, 0.25, na.rm = TRUE),
                q75 = quantile(SPHERICAL_DEPTH, 0.75, na.rm = TRUE),
                .groups = "drop")

    p_ingr_depth_time <- ggplot(ingr_depth_time, aes(x = time_min)) +
      geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#E41A1C", alpha = 0.15) +
      geom_line(aes(y = median_depth), color = "#E41A1C", linewidth = 0.7) +
      labs(title = "Ingressing Cell Depth Over Time",
           subtitle = "Median +/- IQR (absolute time)",
           x = "Time (min)", y = "Spherical depth (um)") +
      theme_pipe()

    p11 <- (p_margin_count | p_margin_frac) /
           (p_ingr_roi_frac | p_roi_vel) /
           (p_roi_depth_comp | p_roi_theta) /
           (p_ingr_heat | p_ingr_depth_time) +
      plot_annotation(
        title = sprintf("%s -- 11 Margin & ROI Density Analysis", sp$name),
        subtitle = "Cell density and behaviour in the ingression zone",
        theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5, color = "grey50")))

    save_pdf(p11, sprintf("%s_11_margin_density.pdf", prefix), width = 16, height = 20)

  }  # end ingression-specific sections

  cat(sprintf("  %s per-species analysis complete.\n", sp$name))

}  # end per-species loop


# #############################################################################
#
#         CROSS-SPECIES COMPARISON DASHBOARD
#
# #############################################################################

cat("\n=== SECTION 12: CROSS-SPECIES COMPARISON ===\n")

if (length(DATA) >= 2) {

  all_tracks <- bind_rows(lapply(DATA, function(d) d$tracks))

  # --- 12a: Normalised nuclei count ---
  norm_counts <- bind_rows(lapply(DATA, function(d) {
    d$spots %>%
      group_by(FRAME) %>%
      summarise(n = n(), .groups = "drop") %>%
      mutate(species  = d$config$name,
             time_min = frame_to_min(FRAME))
  })) %>%
    group_by(species) %>%
    mutate(norm_n    = n / max(n),
           norm_time = (time_min - min(time_min)) / (max(time_min) - min(time_min))) %>%
    ungroup()

  p12a <- ggplot(norm_counts, aes(x = norm_time, y = norm_n, color = species)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = species_colors) +
    labs(title = "Normalised Nuclei Count",
         subtitle = "Both axes [0, 1]",
         x = "Normalised time", y = "Normalised count", color = "Species") +
    theme_pipe()

  # --- 12b: Normalised depth (depth / R) ---
  p12b <- bind_rows(lapply(DATA, function(d) {
    d$spots %>%
      mutate(depth_norm = SPHERICAL_DEPTH / as.numeric(d$sphere["radius"])) %>%
      select(depth_norm, species)
  })) %>%
    ggplot(aes(x = depth_norm, fill = species, color = species)) +
    geom_density(alpha = 0.3, linewidth = 0.6) +
    scale_fill_manual(values = species_colors) +
    scale_color_manual(values = species_colors) +
    labs(title = "Normalised Depth (depth / R)",
         x = "Depth / Radius", y = "Density") +
    theme_pipe()

  # --- 12c: Normalised radial velocity ---
  p12c <- bind_rows(lapply(DATA, function(d) {
    d$spots %>%
      filter(!is.na(RADIAL_VELOCITY_SMOOTH)) %>%
      mutate(vel_norm = RADIAL_VELOCITY_SMOOTH / as.numeric(d$sphere["radius"]) * 1000) %>%
      select(vel_norm, species)
  })) %>%
    ggplot(aes(x = vel_norm, fill = species, color = species)) +
    geom_density(alpha = 0.3, linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = species_colors) +
    scale_color_manual(values = species_colors) +
    coord_cartesian(xlim = c(-10, 10)) +
    labs(title = "Normalised Radial Velocity",
         subtitle = "V_r / R x 1000",
         x = "Normalised velocity (x10^-3 R/frame)", y = "Density") +
    theme_pipe()

  # --- 12d: Instantaneous speed comparison ---
  p12d <- bind_rows(lapply(DATA, function(d) {
    d$spots_vel %>%
      filter(!is.na(inst_speed_um_min), inst_speed_um_min < 20) %>%
      select(inst_speed_um_min, species)
  })) %>%
    ggplot(aes(x = inst_speed_um_min, fill = species, color = species)) +
    geom_density(alpha = 0.3, linewidth = 0.6) +
    scale_fill_manual(values = species_colors) +
    scale_color_manual(values = species_colors) +
    labs(title = "Instantaneous 3D Speed Comparison",
         x = "Speed (um/min)", y = "Density") +
    theme_pipe()

  # --- 12e: Track length comparison ---
  max_frames <- sapply(DATA, function(d) max(d$tracks$frame_end))
  p12e <- all_tracks %>%
    mutate(norm_len = n_spots / max_frames[ifelse(species == "Zebrafish",
                                                   "zebrafish", "medaka")]) %>%
    ggplot(aes(x = norm_len, fill = species, color = species)) +
    geom_density(alpha = 0.3, linewidth = 0.6) +
    scale_fill_manual(values = species_colors) +
    scale_color_manual(values = species_colors) +
    scale_x_log10() +
    labs(title = "Normalised Track Length",
         x = "Track length / total frames (log)", y = "Density") +
    theme_pipe()

  # --- 12f: Depth-stratified velocity comparison ---
  depth_vel_comp <- bind_rows(lapply(DATA, function(d) {
    add_depth_bin(d$spots) %>%
      filter(!is.na(RADIAL_VELOCITY_SMOOTH)) %>%
      group_by(depth_bin) %>%
      summarise(median_vr = median(RADIAL_VELOCITY_SMOOTH),
                n = n(), .groups = "drop") %>%
      mutate(species = d$config$name)
  }))

  p12f <- depth_vel_comp %>%
    ggplot(aes(x = depth_bin, y = median_vr, fill = species)) +
    geom_col(position = "dodge", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    scale_fill_manual(values = species_colors) +
    labs(title = "Median Radial Velocity by Depth",
         x = "Depth layer", y = "Median V_r (um/frame)", fill = "Species") +
    theme_pipe() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # --- 12g: Inward fraction by depth ---
  inward_depth <- bind_rows(lapply(DATA, function(d) {
    add_depth_bin(d$spots) %>%
      filter(!is.na(RADIAL_VELOCITY_SMOOTH)) %>%
      group_by(depth_bin) %>%
      summarise(inward_frac = mean(RADIAL_VELOCITY_SMOOTH < 0),
                n = n(), .groups = "drop") %>%
      mutate(species = d$config$name)
  }))

  p12g <- inward_depth %>%
    ggplot(aes(x = depth_bin, y = inward_frac, fill = species)) +
    geom_col(position = "dodge", alpha = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
    scale_fill_manual(values = species_colors) +
    scale_y_continuous(labels = percent, limits = c(0, 1)) +
    labs(title = "Inward Fraction by Depth",
         x = "Depth layer", y = "Fraction inward", fill = "Species") +
    theme_pipe() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # --- 12h: Turning angle comparison ---
  p12h <- bind_rows(lapply(DATA, function(d) {
    d$spots_vel %>%
      filter(!is.na(turning_angle)) %>%
      select(turning_angle, species)
  })) %>%
    ggplot(aes(x = turning_angle, fill = species, color = species)) +
    geom_density(alpha = 0.3, linewidth = 0.6) +
    geom_vline(xintercept = 90, linetype = "dashed", color = "grey30") +
    scale_fill_manual(values = species_colors) +
    scale_color_manual(values = species_colors) +
    labs(title = "Turning Angle Comparison",
         subtitle = "90 deg = random",
         x = "Turning angle (deg)", y = "Density") +
    theme_pipe()

  p12 <- (p12a | p12b) / (p12c | p12d) / (p12e | p12f) / (p12g | p12h) +
    plot_annotation(
      title = "12 -- Cross-Species Comparison Dashboard",
      subtitle = "Normalised metrics for scale-independent comparison",
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5, color = "grey50")))

  save_pdf(p12, "12_cross_species.pdf", width = 16, height = 20)

} else {
  cat("  Only one species loaded -- skipping cross-species dashboard\n")
}


# #############################################################################
# DONE
# #############################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("  ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n")
cat(sprintf("  Output directory: %s/\n", OUTPUT_DIR))
cat("  Files generated:\n")
for (f in sort(list.files(OUTPUT_DIR, pattern = "\\.pdf$|\\.csv$"))) {
  cat(sprintf("    %s\n", f))
}
cat("\n")
