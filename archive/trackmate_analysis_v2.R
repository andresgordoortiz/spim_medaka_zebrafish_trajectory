# =============================================================================
# TrackMate Analysis — Comparative 4D Embryo Nuclei Tracking Analysis
# =============================================================================
#
# PURPOSE:
#   This is Step 3 of the production pipeline:
#     Step 0: TrackMate → export ALL tracks (no filters)
#     Step 1: embryo_viewer.py → orient, sphere-fit, ROI, ingression detection
#     Step 2: trackmate_filter_and_validate.R → filter + validate (medaka only)
#   → Step 3: THIS SCRIPT → comprehensive biological analysis (both species)
#
# WHAT IT DOES:
#   Loads the enriched CSVs from the napari viewer (Step 1) for both zebrafish
#   and medaka, then performs quality control, depth-stratified analysis,
#   velocity profiling, spatial analysis, and ingression characterization.
#   All plots are structured side-by-side for cross-species comparison.
#
# KEY DESIGN DECISIONS:
#   • NO upfront filtering — all tracks are analysed
#   • Uses SPHERICAL_DEPTH (not Cartesian Z) for depth — more biologically
#     meaningful for a curved embryo surface
#   • Ingression analysis is ZEBRAFISH ONLY (medaka labels are not validated)
#   • All figures have paired panels for future cross-species comparison
#
# INPUTS (per species, from analysis_output_{species}/):
#   - oriented_spots.csv      All spots with enriched spherical coordinates
#   - track_summary.csv       Per-track aggregates
#   - sphere_params.csv       Fitted sphere centre & radius
#   - roi_bounds.csv          ROI bounding box
#   - ingression_params.csv   Ingression detection parameters
#   - analysis_metadata.json  Viewer session metadata
#
# OUTPUTS (in analysis_output/):
#   - 00_data_overview.pdf
#   - 01_nuclei_count_qc.pdf
#   - 02_nuclear_density_voxels.pdf
#   - 03_track_topology.pdf
#   - 04_depth_stratified.pdf
#   - 05_velocity_analysis.pdf
#   - 06_spatial_distribution.pdf
#   - 07_directionality.pdf
#   - 08_intensity_qc.pdf
#   - 09_ingression_analysis.pdf        (zebrafish only)
#   - 10_cross_species_dashboard.pdf
#   - analysis_summary.csv
#
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

# Species configuration
SPECIES <- list(
  zebrafish = list(
    name      = "Zebrafish",
    short     = "ZF",
    color     = "#2166AC",
    data_dir  = "analysis_output_zebrafish",
    has_ingression = TRUE   # validated ingression labels
  ),
  medaka = list(
    name      = "Medaka",
    short     = "MK",
    color     = "#B2182B",
    data_dir  = "analysis_output_medaka",
    filter_dir = "analysis_output",  # filtered data from trackmate_filter_and_validate.R
    has_ingression = FALSE  # ingression labels NOT validated — do not analyse
  )
)

# 3D voxel parameters for density analysis
VOXEL_SIZE_UM <- 50   # µm cube side length

# Spherical depth bins (µm from surface)
DEPTH_BREAKS   <- c(-Inf, 0, 10, 20, 40, 60, Inf)
DEPTH_LABELS   <- c("Outside (<0)", "Surface (0–10)", "Shallow (10–20)",
                     "Mid (20–40)", "Deep (40–60)", "Very deep (>60)")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# =============================================================================
# THEME
# =============================================================================

species_colors <- c("Zebrafish" = "#2166AC", "Medaka" = "#B2182B")

theme_pipe <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1, color = "grey50"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

# Utility: save a ggplot/patchwork to PDF and report
save_pdf <- function(plot, filename, width = 12, height = 8) {
  path <- file.path(OUTPUT_DIR, filename)
  ggsave(path, plot, width = width, height = height, device = cairo_pdf)
  cat(sprintf("  → Saved %s\n", path))
}

# =============================================================================
# HELPER: convert frame to real time
# =============================================================================

frame_to_min  <- function(f) f * FRAME_INTERVAL_MIN
frame_to_hpf  <- function(f) f * FRAME_INTERVAL_SEC / 3600  # hours post-fertilisation offset


cat("\n")
cat(strrep("=", 70), "\n")
cat("  TRACKMATE COMPARATIVE ANALYSIS\n")
cat(strrep("=", 70), "\n")
cat(sprintf("  Frame interval: %d sec | Speed conversion: ×%.0f\n",
            FRAME_INTERVAL_SEC, SPEED_CONVERSION))
cat(sprintf("  Output: %s/\n", OUTPUT_DIR))

# #############################################################################
# SECTION 00: DATA LOADING
# #############################################################################

cat("\n=== SECTION 00: LOADING DATA ===\n")

load_species_data <- function(sp) {
  dir <- sp$data_dir
  if (!dir.exists(dir)) {
    warning(sprintf("Data directory '%s' not found — skipping %s", dir, sp$name))
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
  setDF(spots)  # convert to data.frame for dplyr compatibility
  spots$species <- sp$name
  cat(sprintf("    Spots: %s rows, %d columns, frames %d–%d\n",
              format(nrow(spots), big.mark = ","),
              ncol(spots),
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
  cat(sprintf("    Tracks: %s\n", format(nrow(tsummary), big.mark = ",")))

  # --- Sphere parameters ---
  sphere <- read_csv(file.path(dir, "sphere_params.csv"), show_col_types = FALSE)
  sphere_list <- setNames(sphere$value, sphere$parameter)

  # --- ROI bounds ---
  roi <- read_csv(file.path(dir, "roi_bounds.csv"), show_col_types = FALSE)
  roi_list <- setNames(roi$value, roi$parameter)

  # --- Ingression parameters ---
  ingr <- read_csv(file.path(dir, "ingression_params.csv"), show_col_types = FALSE)
  ingr_list <- setNames(ingr$value, ingr$parameter)

  # --- Metadata ---
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

# Convenience: combined track summaries for cross-species plots
all_tracks <- bind_rows(lapply(DATA, function(d) d$tracks))


# #############################################################################
# SECTION 00b: DATA OVERVIEW TABLE
# #############################################################################

cat("\n=== SECTION 00b: DATA OVERVIEW ===\n")

overview_rows <- lapply(DATA, function(d) {
  sp <- d$spots
  ts <- d$tracks
  tibble(
    Species         = d$config$name,
    `Total spots`   = nrow(sp),
    `Total tracks`  = nrow(ts),
    `Frames`        = length(unique(sp$FRAME)),
    `Frame range`   = sprintf("%d – %d", min(sp$FRAME), max(sp$FRAME)),
    `Duration (min)` = round(max(sp$FRAME) * FRAME_INTERVAL_MIN, 1),
    `Sphere R (µm)` = round(as.numeric(d$sphere["radius"]), 1),
    `Median track length` = median(ts$n_spots),
    `Mean spots/frame` = round(nrow(sp) / length(unique(sp$FRAME)), 0),
    `Depth range (µm)` = sprintf("%.1f – %.1f",
                                  min(sp$SPHERICAL_DEPTH, na.rm = TRUE),
                                  max(sp$SPHERICAL_DEPTH, na.rm = TRUE)),
    `Ingression validated` = d$config$has_ingression,
    `N ingressing tracks` = if (d$config$has_ingression) sum(ts$INGRESSING, na.rm = TRUE) else NA_integer_
  )
})
overview_df <- bind_rows(overview_rows)
print(as_tibble(overview_df))

write_csv(overview_df, file.path(OUTPUT_DIR, "analysis_summary.csv"))
cat("  → Saved analysis_summary.csv\n")

# --- Visual overview ---
p_overview <- overview_df %>%
  select(Species, `Total spots`, `Total tracks`, Frames, `Sphere R (µm)`,
         `Median track length`, `Mean spots/frame`) %>%
  pivot_longer(-Species, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Metric, y = Value, fill = Species)) +
  geom_col(position = "dodge", width = 0.6) +
  scale_fill_manual(values = species_colors) +
  facet_wrap(~Metric, scales = "free", nrow = 2) +
  labs(title = "Dataset Overview", subtitle = "Key metrics per species") +
  theme_pipe() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

save_pdf(p_overview, "00_data_overview.pdf", width = 14, height = 7)


# #############################################################################
# SECTION 01: QC — NUCLEI COUNT OVER TIME
# #############################################################################

cat("\n=== SECTION 01: NUCLEI COUNT OVER TIME ===\n")

# Count spots per frame per species
spots_per_frame <- bind_rows(lapply(DATA, function(d) {
  d$spots %>%
    group_by(FRAME) %>%
    summarise(n_nuclei = n(), .groups = "drop") %>%
    mutate(species  = d$config$name,
           time_min = frame_to_min(FRAME))
}))

p_count <- ggplot(spots_per_frame,
                  aes(x = time_min, y = n_nuclei, color = species)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = species_colors) +
  labs(title = "Detected Nuclei Over Time",
       subtitle = "Total spot count per frame",
       x = "Time (min)", y = "Number of nuclei", color = "Species") +
  theme_pipe()

# Per-species faceted with trend line
p_count_facet <- spots_per_frame %>%
  ggplot(aes(x = time_min, y = n_nuclei)) +
  geom_line(aes(color = species), linewidth = 0.5, alpha = 0.7) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, color = "black",
              linewidth = 0.8) +
  scale_color_manual(values = species_colors) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Nuclei Count with Trend",
       subtitle = "LOESS smoothing (span = 0.3)",
       x = "Time (min)", y = "Number of nuclei") +
  theme_pipe() + theme(legend.position = "none")

# Rate of change
spots_rate <- spots_per_frame %>%
  group_by(species) %>%
  arrange(FRAME) %>%
  mutate(delta = n_nuclei - lag(n_nuclei),
         pct_change = delta / lag(n_nuclei) * 100) %>%
  ungroup()

p_rate <- spots_rate %>%
  filter(!is.na(delta)) %>%
  ggplot(aes(x = time_min, y = delta, color = species)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(alpha = 0.4, linewidth = 0.3) +
  geom_smooth(method = "loess", span = 0.2, se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = species_colors) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Frame-to-Frame Change in Nuclei Count",
       subtitle = "Positive = gain, negative = loss (tracking dropouts or divisions)",
       x = "Time (min)", y = "∆ Nuclei / frame") +
  theme_pipe() + theme(legend.position = "none")

p01 <- p_count / p_count_facet / p_rate +
  plot_annotation(title = "01 — Quality Control: Nuclei Count",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

save_pdf(p01, "01_nuclei_count_qc.pdf", width = 14, height = 14)


# #############################################################################
# SECTION 02: QC — NUCLEAR DENSITY (3D VOXELS) OVER TIME
# #############################################################################

cat("\n=== SECTION 02: NUCLEAR DENSITY (3D VOXELS) ===\n")

compute_voxel_density <- function(spots_df, voxel_size) {
  # Assign each spot to a 3D voxel
  spots_df %>%
    mutate(
      vx = floor(POSITION_X / voxel_size),
      vy = floor(POSITION_Y / voxel_size),
      vz = floor(POSITION_Z / voxel_size)
    ) %>%
    group_by(FRAME, vx, vy, vz) %>%
    summarise(count = n(), .groups = "drop")
}

# Compute per species (may be slow for medaka — use data.table under the hood)
density_stats <- bind_rows(lapply(DATA, function(d) {
  cat(sprintf("    Computing voxel density for %s ...\n", d$config$name))

  dt <- as.data.table(d$spots)[, .(FRAME, POSITION_X, POSITION_Y, POSITION_Z)]
  dt[, `:=`(vx = floor(POSITION_X / VOXEL_SIZE_UM),
            vy = floor(POSITION_Y / VOXEL_SIZE_UM),
            vz = floor(POSITION_Z / VOXEL_SIZE_UM))]

  # Count per voxel per frame
  voxel_counts <- dt[, .(count = .N), by = .(FRAME, vx, vy, vz)]

  # Number of occupied voxels per frame
  n_voxels <- voxel_counts[, .(n_occupied_voxels = .N), by = FRAME]

  # Density statistics per frame
  stats <- voxel_counts[, .(
    mean_density   = mean(count),
    median_density = as.double(median(count)),
    sd_density     = sd(count),
    max_density    = max(count),
    q25_density    = quantile(count, 0.25),
    q75_density    = quantile(count, 0.75),
    total_nuclei   = sum(count)
  ), by = FRAME]

  stats <- merge(stats, n_voxels, by = "FRAME")
  stats$species <- d$config$name
  stats$time_min <- frame_to_min(stats$FRAME)
  as.data.frame(stats)
}))

# Plot: mean density per voxel over time
p_dens_mean <- ggplot(density_stats,
                      aes(x = time_min, y = mean_density, color = species)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = species_colors) +
  labs(title = sprintf("Mean Nuclear Density per Voxel (%dµm³)", VOXEL_SIZE_UM),
       x = "Time (min)", y = "Mean nuclei / voxel", color = "Species") +
  theme_pipe()

# Plot: number of occupied voxels over time (spatial spread)
p_dens_nvox <- density_stats %>%
  ggplot(aes(x = time_min, y = n_occupied_voxels, color = species)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = species_colors) +
  labs(title = "Spatial Spread: Occupied Voxels Over Time",
       subtitle = "Number of non-empty voxels per frame",
       x = "Time (min)", y = "Occupied voxels", color = "Species") +
  theme_pipe()

# Plot: density variability (IQR band)
p_dens_iqr <- density_stats %>%
  ggplot(aes(x = time_min)) +
  geom_ribbon(aes(ymin = q25_density, ymax = q75_density, fill = species),
              alpha = 0.3) +
  geom_line(aes(y = median_density, color = species), linewidth = 0.7) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Voxel Density Distribution Over Time",
       subtitle = "Median (line) ± IQR (ribbon)",
       x = "Time (min)", y = "Nuclei per voxel") +
  theme_pipe() + theme(legend.position = "none")

# Plot: max density (hotspot tracking)
p_dens_max <- density_stats %>%
  ggplot(aes(x = time_min, y = max_density, color = species)) +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = species_colors) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Maximum Voxel Density (Hotspot) Over Time",
       x = "Time (min)", y = "Max nuclei in single voxel") +
  theme_pipe() + theme(legend.position = "none")

p02 <- (p_dens_mean | p_dens_nvox) / p_dens_iqr / p_dens_max +
  plot_annotation(title = "02 — Quality Control: 3D Nuclear Density",
                  subtitle = sprintf("Voxel size: %dµm × %dµm × %dµm",
                                     VOXEL_SIZE_UM, VOXEL_SIZE_UM, VOXEL_SIZE_UM),
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                                plot.subtitle = element_text(hjust = 0.5, color = "grey50")))

save_pdf(p02, "02_nuclear_density_voxels.pdf", width = 14, height = 16)


# #############################################################################
# SECTION 03: QC — TRACK TOPOLOGY
# #############################################################################

cat("\n=== SECTION 03: TRACK TOPOLOGY ===\n")

# Track length distribution
p_tracklen <- all_tracks %>%
  ggplot(aes(x = n_spots, fill = species)) +
  geom_histogram(bins = 80, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = species_colors) +
  scale_x_log10(labels = comma) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Track Length Distribution",
       subtitle = "Number of spots per track (log scale)",
       x = "Track length (spots)", y = "Count") +
  theme_pipe() + theme(legend.position = "none")

# Track temporal coverage
p_coverage <- all_tracks %>%
  mutate(duration = frame_to_min(frame_end - frame_start),
         coverage_pct = n_spots / (frame_end - frame_start + 1) * 100) %>%
  ggplot(aes(x = duration, y = coverage_pct, color = species)) +
  geom_point(alpha = 0.02, size = 0.3) +
  geom_density2d(linewidth = 0.4) +
  scale_color_manual(values = species_colors) +
  facet_wrap(~species, ncol = 2) +
  labs(title = "Track Duration vs. Temporal Coverage",
       subtitle = "Coverage = spots / expected frames × 100%",
       x = "Duration (min)", y = "Coverage (%)") +
  coord_cartesian(ylim = c(0, 105)) +
  theme_pipe() + theme(legend.position = "none")

# Track start/end frame distribution
track_timing <- all_tracks %>%
  select(species, frame_start, frame_end) %>%
  pivot_longer(cols = c(frame_start, frame_end),
               names_to = "event", values_to = "frame") %>%
  mutate(event = ifelse(event == "frame_start", "Track starts", "Track ends"),
         time_min = frame_to_min(frame))

p_timing <- track_timing %>%
  ggplot(aes(x = time_min, fill = event)) +
  geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("Track starts" = "#66C2A5", "Track ends" = "#FC8D62")) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Track Birth / Death Over Time",
       subtitle = "When tracks begin and end",
       x = "Time (min)", y = "Count") +
  theme_pipe()

# Active tracks over time (tracks alive at each frame)
active_per_frame <- bind_rows(lapply(DATA, function(d) {
  frames <- seq(min(d$tracks$frame_start), max(d$tracks$frame_end))
  tibble(
    FRAME = frames,
    n_active = sapply(frames, function(f) {
      sum(d$tracks$frame_start <= f & d$tracks$frame_end >= f)
    }),
    species = d$config$name,
    time_min = frame_to_min(frames)
  )
}))

p_active <- active_per_frame %>%
  ggplot(aes(x = time_min, y = n_active, color = species)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = species_colors) +
  labs(title = "Active Tracks Over Time",
       subtitle = "Tracks alive at each frame",
       x = "Time (min)", y = "Active tracks", color = "Species") +
  theme_pipe()

p03 <- (p_tracklen / p_coverage) | (p_timing / p_active)

p03 <- p03 +
  plot_annotation(title = "03 — Track Topology",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

save_pdf(p03, "03_track_topology.pdf", width = 16, height = 12)


# #############################################################################
# SECTION 04: DEPTH-STRATIFIED ANALYSIS (SPHERICAL DEPTH)
# #############################################################################

cat("\n=== SECTION 04: DEPTH-STRATIFIED ANALYSIS ===\n")

# Assign depth bins to all spots
add_depth_bin <- function(spots_df) {
  spots_df %>%
    mutate(depth_bin = cut(SPHERICAL_DEPTH,
                           breaks = DEPTH_BREAKS,
                           labels = DEPTH_LABELS,
                           ordered_result = TRUE))
}

# Nuclei per depth bin per frame per species
depth_counts <- bind_rows(lapply(DATA, function(d) {
  add_depth_bin(d$spots) %>%
    group_by(FRAME, depth_bin, species) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(time_min = frame_to_min(FRAME))
}))

# Stacked area: nuclei distribution across depth over time
p_depth_area <- depth_counts %>%
  ggplot(aes(x = time_min, y = n, fill = depth_bin)) +
  geom_area(position = "stack", alpha = 0.8) +
  scale_fill_viridis_d(option = "mako", direction = -1) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Nuclei Distribution Across Depth Layers Over Time",
       subtitle = "Spherical depth from fitted surface (µm)",
       x = "Time (min)", y = "Number of nuclei", fill = "Depth layer") +
  theme_pipe()

# Proportion version
depth_props <- depth_counts %>%
  group_by(FRAME, species) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_depth_prop <- depth_props %>%
  ggplot(aes(x = time_min, y = prop, fill = depth_bin)) +
  geom_area(position = "stack", alpha = 0.8) +
  scale_fill_viridis_d(option = "mako", direction = -1) +
  scale_y_continuous(labels = percent) +
  facet_wrap(~species, ncol = 2) +
  labs(title = "Proportional Depth Distribution Over Time",
       x = "Time (min)", y = "Fraction of nuclei", fill = "Depth layer") +
  theme_pipe()

# Depth distribution comparison: violin per species
p_depth_violin <- bind_rows(lapply(DATA, function(d) {
  d$spots %>% select(SPHERICAL_DEPTH, species)
})) %>%
  ggplot(aes(x = species, y = SPHERICAL_DEPTH, fill = species)) +
  geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(values = species_colors) +
  labs(title = "Overall Depth Distribution",
       x = "", y = "Spherical depth (µm)") +
  theme_pipe() + theme(legend.position = "none")

# Depth over time: heatmap (2D histogram)
p_depth_heat <- bind_rows(lapply(DATA, function(d) {
  d$spots %>%
    select(FRAME, SPHERICAL_DEPTH, species) %>%
    mutate(time_min = frame_to_min(FRAME))
})) %>%
  ggplot(aes(x = time_min, y = SPHERICAL_DEPTH)) +
  geom_bin2d(bins = c(80, 60)) +
  scale_fill_viridis_c(option = "inferno", trans = "log10",
                       labels = comma) +
  facet_wrap(~species, scales = "free_x", ncol = 2) +
  labs(title = "Depth × Time Heatmap",
       subtitle = "Log₁₀ count of nuclei",
       x = "Time (min)", y = "Spherical depth (µm)", fill = "Count") +
  theme_pipe()

p04 <- p_depth_area / p_depth_prop / (p_depth_violin | p_depth_heat) +
  plot_annotation(title = "04 — Depth-Stratified Analysis",
                  subtitle = "Using spherical depth from fitted embryo surface",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                                plot.subtitle = element_text(hjust = 0.5, color = "grey50")))

save_pdf(p04, "04_depth_stratified.pdf", width = 16, height = 18)


# #############################################################################
# SECTION 05: VELOCITY ANALYSIS
# #############################################################################

cat("\n=== SECTION 05: VELOCITY ANALYSIS ===\n")

# Radial velocity distributions
p_vel_hist <- bind_rows(lapply(DATA, function(d) {
  d$spots %>%
    filter(!is.na(RADIAL_VELOCITY_SMOOTH)) %>%
    select(RADIAL_VELOCITY_SMOOTH, species)
})) %>%
  ggplot(aes(x = RADIAL_VELOCITY_SMOOTH, fill = species)) +
  geom_histogram(bins = 100, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  scale_fill_manual(values = species_colors) +
  coord_cartesian(xlim = c(-5, 5)) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Radial Velocity Distribution",
       subtitle = "Negative = inward (toward interior), Positive = outward",
       x = "Smoothed radial velocity (µm/frame)", y = "Count") +
  theme_pipe() + theme(legend.position = "none")

# Mean radial velocity over time
vel_over_time <- bind_rows(lapply(DATA, function(d) {
  d$spots %>%
    filter(!is.na(RADIAL_VELOCITY_SMOOTH)) %>%
    group_by(FRAME) %>%
    summarise(
      mean_vr   = mean(RADIAL_VELOCITY_SMOOTH),
      median_vr = median(RADIAL_VELOCITY_SMOOTH),
      q25_vr    = quantile(RADIAL_VELOCITY_SMOOTH, 0.25),
      q75_vr    = quantile(RADIAL_VELOCITY_SMOOTH, 0.75),
      .groups   = "drop"
    ) %>%
    mutate(species = d$config$name,
           time_min = frame_to_min(FRAME))
}))

p_vel_time <- vel_over_time %>%
  ggplot(aes(x = time_min)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = q25_vr, ymax = q75_vr, fill = species), alpha = 0.2) +
  geom_line(aes(y = median_vr, color = species), linewidth = 0.7) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  facet_wrap(~species, scales = "free", ncol = 2) +
  labs(title = "Radial Velocity Over Time",
       subtitle = "Median (line) ± IQR (ribbon)",
       x = "Time (min)", y = "Radial velocity (µm/frame)") +
  theme_pipe() + theme(legend.position = "none")

# Velocity by depth layer
vel_by_depth <- bind_rows(lapply(DATA, function(d) {
  add_depth_bin(d$spots) %>%
    filter(!is.na(RADIAL_VELOCITY_SMOOTH)) %>%
    select(RADIAL_VELOCITY_SMOOTH, depth_bin, species)
}))

p_vel_depth <- vel_by_depth %>%
  ggplot(aes(x = depth_bin, y = RADIAL_VELOCITY_SMOOTH, fill = species)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.7, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_fill_manual(values = species_colors) +
  coord_cartesian(ylim = c(-3, 3)) +
  facet_wrap(~species, ncol = 2) +
  labs(title = "Radial Velocity by Depth Layer",
       x = "Depth layer", y = "Radial velocity (µm/frame)") +
  theme_pipe() + theme(legend.position = "none",
                       axis.text.x = element_text(angle = 30, hjust = 1))

# Per-track velocity: median V_r distribution
p_track_vel <- all_tracks %>%
  filter(!is.na(TRACK_MEDIAN_RADIAL_VEL)) %>%
  ggplot(aes(x = TRACK_MEDIAN_RADIAL_VEL, fill = species)) +
  geom_histogram(bins = 80, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  scale_fill_manual(values = species_colors) +
  coord_cartesian(xlim = c(-3, 3)) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Per-Track Median Radial Velocity",
       x = "Median radial velocity (µm/frame)", y = "Number of tracks") +
  theme_pipe() + theme(legend.position = "none")

p05 <- p_vel_hist / p_vel_time / p_vel_depth / p_track_vel +
  plot_annotation(title = "05 — Velocity Analysis",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

save_pdf(p05, "05_velocity_analysis.pdf", width = 16, height = 20)


# #############################################################################
# SECTION 06: SPATIAL DISTRIBUTION & ANGULAR COVERAGE
# #############################################################################

cat("\n=== SECTION 06: SPATIAL DISTRIBUTION ===\n")

# Angular distribution: theta vs phi (spherical coordinates)
p_angular <- bind_rows(lapply(DATA, function(d) {
  d$spots %>%
    select(THETA_DEG, PHI_DEG, species) %>%
    filter(!is.na(THETA_DEG), !is.na(PHI_DEG))
})) %>%
  ggplot(aes(x = PHI_DEG, y = THETA_DEG)) +
  geom_bin2d(bins = c(60, 40)) +
  scale_fill_viridis_c(option = "plasma", trans = "log10", labels = comma) +
  facet_wrap(~species, ncol = 2) +
  labs(title = "Angular Coverage (Spherical Coordinates)",
       subtitle = "Θ (polar angle) vs Φ (azimuthal angle)",
       x = "Φ (azimuthal, deg)", y = "Θ (polar, deg)", fill = "Count") +
  theme_pipe()

# XY projection of all spots (2D density)
p_xy <- bind_rows(lapply(DATA, function(d) {
  d$spots %>% select(POSITION_X, POSITION_Y, species)
})) %>%
  ggplot(aes(x = POSITION_X, y = POSITION_Y)) +
  geom_bin2d(bins = 80) +
  scale_fill_viridis_c(option = "plasma", trans = "log10", labels = comma) +
  facet_wrap(~species, scales = "free", ncol = 2) +
  labs(title = "XY Projection of All Nuclei",
       x = "X (µm)", y = "Y (µm)", fill = "Count") +
  theme_pipe()

# Track centroid positions
p_centroids <- all_tracks %>%
  ggplot(aes(x = mean_x, y = mean_y, color = species)) +
  geom_point(alpha = 0.05, size = 0.3) +
  scale_color_manual(values = species_colors) +
  facet_wrap(~species, scales = "free", ncol = 2) +
  labs(title = "Track Centroid Positions (XY)",
       x = "Mean X (µm)", y = "Mean Y (µm)") +
  theme_pipe() + theme(legend.position = "none")

# Depth vs angular position
p_depth_angle <- bind_rows(lapply(DATA, function(d) {
  d$spots %>%
    select(THETA_DEG, SPHERICAL_DEPTH, species) %>%
    filter(!is.na(THETA_DEG))
})) %>%
  ggplot(aes(x = THETA_DEG, y = SPHERICAL_DEPTH)) +
  geom_bin2d(bins = c(50, 50)) +
  scale_fill_viridis_c(option = "inferno", trans = "log10", labels = comma) +
  facet_wrap(~species, ncol = 2) +
  labs(title = "Depth vs. Polar Angle",
       x = "Θ (polar, deg)", y = "Spherical depth (µm)", fill = "Count") +
  theme_pipe()

p06 <- (p_angular / p_xy) | (p_centroids / p_depth_angle)

p06 <- p06 +
  plot_annotation(title = "06 — Spatial Distribution & Angular Coverage",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

save_pdf(p06, "06_spatial_distribution.pdf", width = 16, height = 14)


# #############################################################################
# SECTION 07: DIRECTIONALITY & MOVEMENT PATTERNS
# #############################################################################

cat("\n=== SECTION 07: DIRECTIONALITY ===\n")

# Inward fraction distribution per track
p_inward <- all_tracks %>%
  filter(!is.na(TRACK_INWARD_FRACTION)) %>%
  ggplot(aes(x = TRACK_INWARD_FRACTION, fill = species)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey30") +
  scale_fill_manual(values = species_colors) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Inward Fraction per Track",
       subtitle = ">0.5 = predominantly inward movement",
       x = "Fraction of frames with inward velocity", y = "Count") +
  theme_pipe() + theme(legend.position = "none")

# Net radial displacement
p_netdisp <- all_tracks %>%
  filter(!is.na(TRACK_NET_RADIAL_DISP)) %>%
  ggplot(aes(x = TRACK_NET_RADIAL_DISP, fill = species)) +
  geom_histogram(bins = 80, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  scale_fill_manual(values = species_colors) +
  coord_cartesian(xlim = quantile(all_tracks$TRACK_NET_RADIAL_DISP,
                                  c(0.01, 0.99), na.rm = TRUE)) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Net Radial Displacement per Track",
       subtitle = "Negative = net inward, Positive = net outward",
       x = "Net radial displacement (µm)", y = "Count") +
  theme_pipe() + theme(legend.position = "none")

# Sustained inward movement
p_sustained <- all_tracks %>%
  filter(!is.na(TRACK_SUSTAINED_INWARD), TRACK_SUSTAINED_INWARD > 0) %>%
  ggplot(aes(x = TRACK_SUSTAINED_INWARD, fill = species)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = species_colors) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Max Sustained Consecutive Inward Frames",
       subtitle = "Tracks with ≥1 sustained inward frame",
       x = "Consecutive inward frames", y = "Count") +
  theme_pipe() + theme(legend.position = "none")

# Net displacement vs track length — coloured by depth
p_disp_len <- all_tracks %>%
  filter(!is.na(TRACK_NET_RADIAL_DISP)) %>%
  ggplot(aes(x = n_spots, y = TRACK_NET_RADIAL_DISP, color = mean_SPHERICAL_DEPTH)) +
  geom_point(alpha = 0.08, size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_viridis_c(option = "mako", direction = -1) +
  scale_x_log10() +
  facet_wrap(~species, scales = "free", ncol = 2) +
  labs(title = "Net Radial Displacement vs. Track Length",
       subtitle = "Coloured by mean spherical depth",
       x = "Track length (spots, log)", y = "Net radial displacement (µm)",
       color = "Depth (µm)") +
  theme_pipe()

p07 <- (p_inward / p_netdisp) | (p_sustained / p_disp_len)

p07 <- p07 +
  plot_annotation(title = "07 — Directionality & Movement Patterns",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

save_pdf(p07, "07_directionality.pdf", width = 16, height = 14)


# #############################################################################
# SECTION 08: INTENSITY & SIGNAL QUALITY
# #############################################################################

cat("\n=== SECTION 08: INTENSITY & SIGNAL QUALITY ===\n")

# Intensity over time
intensity_time <- bind_rows(lapply(DATA, function(d) {
  d$spots %>%
    filter(!is.na(MEAN_INTENSITY_CH1)) %>%
    group_by(FRAME) %>%
    summarise(
      mean_int   = mean(MEAN_INTENSITY_CH1),
      median_int = median(MEAN_INTENSITY_CH1),
      mean_snr   = mean(SNR_CH1, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(species = d$config$name,
           time_min = frame_to_min(FRAME))
}))

p_int_time <- intensity_time %>%
  ggplot(aes(x = time_min, y = median_int, color = species)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = species_colors) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Median Spot Intensity Over Time",
       x = "Time (min)", y = "Median intensity (CH1)") +
  theme_pipe() + theme(legend.position = "none")

p_snr_time <- intensity_time %>%
  ggplot(aes(x = time_min, y = mean_snr, color = species)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = species_colors) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Mean SNR Over Time",
       subtitle = "Signal-to-noise ratio",
       x = "Time (min)", y = "Mean SNR (CH1)") +
  theme_pipe() + theme(legend.position = "none")

# Intensity vs depth
p_int_depth <- bind_rows(lapply(DATA, function(d) {
  d$spots %>%
    filter(!is.na(MEAN_INTENSITY_CH1)) %>%
    select(SPHERICAL_DEPTH, MEAN_INTENSITY_CH1, SNR_CH1, species)
})) %>%
  ggplot(aes(x = SPHERICAL_DEPTH, y = MEAN_INTENSITY_CH1)) +
  geom_bin2d(bins = c(50, 50)) +
  scale_fill_viridis_c(option = "magma", trans = "log10", labels = comma) +
  facet_wrap(~species, scales = "free", ncol = 2) +
  labs(title = "Intensity vs. Spherical Depth",
       subtitle = "Deeper nuclei may have lower signal due to scattering",
       x = "Spherical depth (µm)", y = "Mean intensity (CH1)", fill = "Count") +
  theme_pipe()

# Quality score distribution
p_quality <- bind_rows(lapply(DATA, function(d) {
  d$spots %>% select(QUALITY, species)
})) %>%
  ggplot(aes(x = QUALITY, fill = species)) +
  geom_histogram(bins = 80, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = species_colors) +
  scale_x_log10(labels = comma) +
  facet_wrap(~species, scales = "free_y", ncol = 2) +
  labs(title = "Spot Quality Score Distribution",
       x = "Quality (log scale)", y = "Count") +
  theme_pipe() + theme(legend.position = "none")

p08 <- p_int_time / p_snr_time / p_int_depth / p_quality +
  plot_annotation(title = "08 — Intensity & Signal Quality",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

save_pdf(p08, "08_intensity_qc.pdf", width = 14, height = 18)


# #############################################################################
# SECTION 09: INGRESSION ANALYSIS (ZEBRAFISH ONLY)
# #############################################################################

cat("\n=== SECTION 09: INGRESSION ANALYSIS ===\n")

# Only run for species with validated ingression labels
ingression_species <- names(DATA)[sapply(DATA, function(d) d$config$has_ingression)]

if (length(ingression_species) == 0) {
  cat("  ⚠ No species with validated ingression → skipping section 09\n")
} else {
  for (sp_key in ingression_species) {
    d <- DATA[[sp_key]]
    sp_name <- d$config$name
    cat(sprintf("  Analysing ingression for %s ...\n", sp_name))

    spots  <- d$spots
    tracks <- d$tracks

    ingr_tracks <- tracks %>% filter(INGRESSING == TRUE | INGRESSING == "True")
    non_ingr_tracks <- tracks %>% filter(INGRESSING == FALSE | INGRESSING == "False")

    n_ingr <- nrow(ingr_tracks)
    n_total <- nrow(tracks)
    cat(sprintf("    Ingressing: %d / %d tracks (%.1f%%)\n",
                n_ingr, n_total, 100 * n_ingr / n_total))

    # Mark spots as ingressing
    ingr_ids <- ingr_tracks$TRACK_ID
    spots_ingr <- spots %>%
      mutate(is_ingressing_track = TRACK_ID %in% ingr_ids)

    # --- 9a: Ingression onset timing ---
    p_onset <- ingr_tracks %>%
      filter(!is.na(INGRESSION_ONSET_FRAME)) %>%
      mutate(onset_min = frame_to_min(INGRESSION_ONSET_FRAME)) %>%
      ggplot(aes(x = onset_min)) +
      geom_histogram(bins = 40, fill = "#E41A1C", alpha = 0.7) +
      labs(title = sprintf("Ingression Onset Timing (%s)", sp_name),
           subtitle = sprintf("%d ingressing tracks", n_ingr),
           x = "Onset time (min)", y = "Number of tracks") +
      theme_pipe()

    # --- 9b: Where do ingressing cells start? (depth at onset) ---
    ingr_onset_depth <- spots_ingr %>%
      filter(is_ingressing_track,
             !is.na(INGRESSION_ONSET_FRAME)) %>%
      group_by(TRACK_ID) %>%
      filter(FRAME == first(INGRESSION_ONSET_FRAME)) %>%
      ungroup()

    p_onset_depth <- ingr_onset_depth %>%
      ggplot(aes(x = SPHERICAL_DEPTH)) +
      geom_histogram(bins = 30, fill = "#E41A1C", alpha = 0.7) +
      labs(title = "Depth at Ingression Onset",
           subtitle = "Spherical depth when inward movement begins",
           x = "Spherical depth (µm)", y = "Count") +
      theme_pipe()

    # --- 9c: Spatial distribution of ingressing vs non-ingressing ---
    p_ingr_spatial <- spots_ingr %>%
      ggplot(aes(x = THETA_DEG, y = SPHERICAL_DEPTH,
                 color = is_ingressing_track)) +
      geom_point(data = . %>% filter(!is_ingressing_track),
                 alpha = 0.005, size = 0.1, color = "grey70") +
      geom_point(data = . %>% filter(is_ingressing_track),
                 alpha = 0.05, size = 0.3, color = "#E41A1C") +
      labs(title = "Ingressing vs. Non-Ingressing Nuclei",
           subtitle = "Θ-depth plane (ingressing in red)",
           x = "Θ (polar angle, deg)", y = "Spherical depth (µm)") +
      theme_pipe() + theme(legend.position = "none")

    # --- 9d: Velocity profile around ingression onset ---
    # Align all ingressing tracks to their onset frame
    onset_frames <- ingr_tracks %>%
      filter(!is.na(INGRESSION_ONSET_FRAME)) %>%
      select(TRACK_ID, onset = INGRESSION_ONSET_FRAME)

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
        n         = n(),
        .groups   = "drop"
      ) %>%
      filter(n >= 5)  # require at least 5 tracks contributing

    p_vel_profile <- vel_profile %>%
      filter(abs(rel_frame) <= 30) %>%
      ggplot(aes(x = rel_frame * FRAME_INTERVAL_MIN, y = mean_vr)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#E41A1C",
                 linewidth = 0.8) +
      geom_ribbon(aes(ymin = mean_vr - se_vr, ymax = mean_vr + se_vr),
                  fill = "#E41A1C", alpha = 0.2) +
      geom_line(linewidth = 0.8, color = "#E41A1C") +
      annotate("text", x = 0.5, y = max(vel_profile$mean_vr) * 0.8,
               label = "Onset", color = "#E41A1C", hjust = 0, fontface = "italic") +
      labs(title = "Radial Velocity Aligned to Ingression Onset",
           subtitle = "Mean ± SE across all ingressing tracks",
           x = "Time relative to onset (min)", y = "Radial velocity (µm/frame)") +
      theme_pipe()

    # --- 9e: Depth trajectory of ingressing cells ---
    depth_traj <- spots %>%
      filter(TRACK_ID %in% ingr_ids) %>%
      inner_join(onset_frames, by = "TRACK_ID") %>%
      mutate(rel_frame = FRAME - onset)

    depth_traj_summary <- depth_traj %>%
      group_by(rel_frame) %>%
      summarise(
        mean_depth   = mean(SPHERICAL_DEPTH, na.rm = TRUE),
        median_depth = median(SPHERICAL_DEPTH, na.rm = TRUE),
        q25_depth    = quantile(SPHERICAL_DEPTH, 0.25, na.rm = TRUE),
        q75_depth    = quantile(SPHERICAL_DEPTH, 0.75, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>%
      filter(n >= 5)

    p_depth_traj <- depth_traj_summary %>%
      filter(abs(rel_frame) <= 30) %>%
      ggplot(aes(x = rel_frame * FRAME_INTERVAL_MIN)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#E41A1C",
                 linewidth = 0.8) +
      geom_ribbon(aes(ymin = q25_depth, ymax = q75_depth),
                  fill = "#E41A1C", alpha = 0.15) +
      geom_line(aes(y = median_depth), linewidth = 0.8, color = "#E41A1C") +
      labs(title = "Depth Trajectory Around Ingression Onset",
           subtitle = "Median ± IQR of spherical depth",
           x = "Time relative to onset (min)", y = "Spherical depth (µm)") +
      theme_pipe()

    # --- 9f: Ingression score distribution ---
    p_score <- tracks %>%
      filter(!is.na(INGRESSION_SCORE)) %>%
      ggplot(aes(x = INGRESSION_SCORE, fill = INGRESSING)) +
      geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
      scale_fill_manual(values = c("TRUE" = "#E41A1C", "True" = "#E41A1C",
                                   "FALSE" = "grey60", "False" = "grey60"),
                        labels = c("Non-ingressing", "Ingressing")) +
      labs(title = "Ingression Score Distribution",
           subtitle = "Weighted composite score (0 = none, 1 = strong)",
           x = "Ingression score", y = "Number of tracks", fill = "") +
      theme_pipe()

    # --- 9g: Ingressing track properties comparison ---
    comp_data <- tracks %>%
      mutate(group = ifelse(INGRESSING == TRUE | INGRESSING == "True",
                            "Ingressing", "Non-ingressing"))

    p_comp_len <- comp_data %>%
      ggplot(aes(x = group, y = n_spots, fill = group)) +
      geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
      scale_fill_manual(values = c("Ingressing" = "#E41A1C",
                                   "Non-ingressing" = "grey60")) +
      labs(title = "Track Length", y = "Number of spots", x = "") +
      theme_pipe() + theme(legend.position = "none")

    p_comp_depth <- comp_data %>%
      ggplot(aes(x = group, y = mean_SPHERICAL_DEPTH, fill = group)) +
      geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
      scale_fill_manual(values = c("Ingressing" = "#E41A1C",
                                   "Non-ingressing" = "grey60")) +
      labs(title = "Mean Depth", y = "Spherical depth (µm)", x = "") +
      theme_pipe() + theme(legend.position = "none")

    p_comp_vel <- comp_data %>%
      filter(!is.na(TRACK_MEDIAN_RADIAL_VEL)) %>%
      ggplot(aes(x = group, y = TRACK_MEDIAN_RADIAL_VEL, fill = group)) +
      geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
      scale_fill_manual(values = c("Ingressing" = "#E41A1C",
                                   "Non-ingressing" = "grey60")) +
      coord_cartesian(ylim = c(-3, 3)) +
      labs(title = "Median Radial Velocity", y = "µm/frame", x = "") +
      theme_pipe() + theme(legend.position = "none")

    p_comps <- p_comp_len | p_comp_depth | p_comp_vel

    # --- 9h: Ingression rate over time ---
    # How many tracks start ingressing per time window?
    ingr_rate <- ingr_tracks %>%
      filter(!is.na(INGRESSION_ONSET_FRAME)) %>%
      mutate(onset_min = frame_to_min(INGRESSION_ONSET_FRAME),
             time_bin = floor(onset_min / 5) * 5) %>%  # 5-min bins
      group_by(time_bin) %>%
      summarise(n_new = n(), .groups = "drop")

    p_rate <- ingr_rate %>%
      ggplot(aes(x = time_bin, y = n_new)) +
      geom_col(fill = "#E41A1C", alpha = 0.7) +
      labs(title = "Ingression Rate Over Time",
           subtitle = "New ingression onsets per 5-min window",
           x = "Time (min)", y = "New ingressing tracks") +
      theme_pipe()

    # --- Assemble ---
    p09_top    <- (p_onset | p_onset_depth | p_score)
    p09_mid    <- (p_ingr_spatial | p_vel_profile)
    p09_bot    <- (p_depth_traj | p_rate)

    p09 <- p09_top / p09_mid / p09_bot / p_comps +
      plot_annotation(
        title = sprintf("09 — Ingression Analysis (%s)", sp_name),
        subtitle = sprintf("%d ingressing tracks out of %s total (%.1f%%)",
                           n_ingr, format(n_total, big.mark = ","),
                           100 * n_ingr / n_total),
        theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5, color = "grey50"))
      )

    save_pdf(p09, sprintf("09_ingression_analysis_%s.pdf", tolower(sp_name)),
             width = 18, height = 22)
  }
}


# #############################################################################
# SECTION 10: CROSS-SPECIES COMPARISON DASHBOARD
# #############################################################################

cat("\n=== SECTION 10: CROSS-SPECIES DASHBOARD ===\n")

if (length(DATA) >= 2) {

  # 10a: Normalised nuclei count (fraction of max per species)
  norm_counts <- spots_per_frame %>%
    group_by(species) %>%
    mutate(norm_n = n_nuclei / max(n_nuclei)) %>%
    ungroup()

  # Normalise time axis to [0, 1] for comparison
  norm_counts <- norm_counts %>%
    group_by(species) %>%
    mutate(norm_time = (time_min - min(time_min)) / (max(time_min) - min(time_min))) %>%
    ungroup()

  p10a <- norm_counts %>%
    ggplot(aes(x = norm_time, y = norm_n, color = species)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = species_colors) +
    labs(title = "Normalised Nuclei Count Over Normalised Time",
         subtitle = "Both axes scaled [0, 1] for comparison",
         x = "Normalised time", y = "Normalised nuclei count",
         color = "Species") +
    theme_pipe()

  # 10b: Normalised depth distribution (depth / sphere radius)
  p10b <- bind_rows(lapply(DATA, function(d) {
    d$spots %>%
      mutate(depth_norm = SPHERICAL_DEPTH / as.numeric(d$sphere["radius"])) %>%
      select(depth_norm, species)
  })) %>%
    ggplot(aes(x = depth_norm, fill = species, color = species)) +
    geom_density(alpha = 0.3, linewidth = 0.6) +
    scale_fill_manual(values = species_colors) +
    scale_color_manual(values = species_colors) +
    labs(title = "Normalised Depth Distribution",
         subtitle = "Depth / sphere radius — scale-independent comparison",
         x = "Normalised depth (depth / R)", y = "Density",
         fill = "Species", color = "Species") +
    theme_pipe()

  # 10c: Normalised velocity distribution (velocity / sphere radius)
  p10c <- bind_rows(lapply(DATA, function(d) {
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
    labs(title = "Normalised Radial Velocity Distribution",
         subtitle = "Velocity / R × 1000 — scale-independent comparison",
         x = "Normalised velocity (×10⁻³ R/frame)", y = "Density",
         fill = "Species", color = "Species") +
    theme_pipe()

  # 10d: Track length distribution (normalised to experiment duration)
  max_frames <- sapply(DATA, function(d) max(d$tracks$frame_end))
  p10d <- all_tracks %>%
    mutate(norm_len = n_spots / max_frames[ifelse(species == "Zebrafish",
                                                   "zebrafish", "medaka")]) %>%
    ggplot(aes(x = norm_len, fill = species, color = species)) +
    geom_density(alpha = 0.3, linewidth = 0.6) +
    scale_fill_manual(values = species_colors) +
    scale_color_manual(values = species_colors) +
    scale_x_log10() +
    labs(title = "Normalised Track Length",
         subtitle = "Track length / total frames — coverage fraction",
         x = "Normalised track length (log)", y = "Density") +
    theme_pipe()

  # 10e: Depth-stratified velocity comparison
  depth_vel_comp <- bind_rows(lapply(DATA, function(d) {
    add_depth_bin(d$spots) %>%
      filter(!is.na(RADIAL_VELOCITY_SMOOTH)) %>%
      group_by(depth_bin) %>%
      summarise(
        mean_vr   = mean(RADIAL_VELOCITY_SMOOTH),
        median_vr = median(RADIAL_VELOCITY_SMOOTH),
        sd_vr     = sd(RADIAL_VELOCITY_SMOOTH),
        n         = n(),
        .groups   = "drop"
      ) %>%
      mutate(species = d$config$name)
  }))

  p10e <- depth_vel_comp %>%
    ggplot(aes(x = depth_bin, y = median_vr, fill = species)) +
    geom_col(position = "dodge", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    scale_fill_manual(values = species_colors) +
    labs(title = "Median Radial Velocity by Depth Layer",
         subtitle = "Cross-species comparison of movement direction per layer",
         x = "Depth layer", y = "Median radial velocity (µm/frame)",
         fill = "Species") +
    theme_pipe() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # 10f: Inward fraction by depth
  inward_depth <- bind_rows(lapply(DATA, function(d) {
    add_depth_bin(d$spots) %>%
      filter(!is.na(RADIAL_VELOCITY_SMOOTH)) %>%
      group_by(depth_bin) %>%
      summarise(
        inward_frac = mean(RADIAL_VELOCITY_SMOOTH < 0),
        n = n(),
        .groups = "drop"
      ) %>%
      mutate(species = d$config$name)
  }))

  p10f <- inward_depth %>%
    ggplot(aes(x = depth_bin, y = inward_frac, fill = species)) +
    geom_col(position = "dodge", alpha = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
    scale_fill_manual(values = species_colors) +
    scale_y_continuous(labels = percent, limits = c(0, 1)) +
    labs(title = "Fraction of Inward-Moving Nuclei by Depth",
         x = "Depth layer", y = "Fraction inward", fill = "Species") +
    theme_pipe() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  p10 <- (p10a | p10b) / (p10c | p10d) / (p10e | p10f) +
    plot_annotation(
      title = "10 — Cross-Species Comparison Dashboard",
      subtitle = "All metrics normalised for scale-independent comparison",
      theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5, color = "grey50"))
    )

  save_pdf(p10, "10_cross_species_dashboard.pdf", width = 16, height = 18)

} else {
  cat("  ⚠ Only one species loaded — skipping cross-species dashboard\n")
}


# #############################################################################
# DONE
# #############################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("  ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n")
cat(sprintf("  Output directory: %s/\n", OUTPUT_DIR))
cat(sprintf("  Files generated:\n"))
for (f in sort(list.files(OUTPUT_DIR))) {
  cat(sprintf("    • %s\n", f))
}
cat("\n")
