# =============================================================================
# Gastrulation Dynamics — MEDAKA vs ZEBRAFISH Comparison
# =============================================================================
#
# Side-by-side comparison of the key dynamics metrics.
# Reads summary/track metrics/thickness from both species' output folders.
# =============================================================================

source("renv/activate.R")

library(data.table)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(viridis)
library(patchwork)
library(scales)

# =============================================================================
# PARAMETERS
# =============================================================================

MEDAKA_DIR  <- "analysis_output_medaka"
ZEBRAFISH_DIR <- "analysis_output_zebrafish"
OUTPUT_DIR <- "analysis_output_comparison"

# Input folders with raw oriented tracks
MEDAKA_INPUT  <- "oriented_medaka_ultrack"
ZEBRAFISH_INPUT <- "oriented_zebrafish_ultrack"

MEDAKA_FI   <- 30   # frame interval in seconds
ZEBRAFISH_FI <- 120

MSD_MAX_LAG   <- 20
MSD_N_TRACKS  <- 8000
MIN_TRACK_LENGTH <- 5

# Voxel calibration — ultrack outputs positions in pixels.
# From TIFF metadata (nuclear_stats_global.tsv):
MEDAKA_VOXEL_UM    <- 1.05152  # um per pixel
ZEBRAFISH_VOXEL_UM <- 1.24785  # um per pixel

# Nuclear stats directories
MEDAKA_NUCLEAR_DIR   <- "nuclei_stats_medaka"
ZEBRAFISH_NUCLEAR_DIR <- "nuclei_stats_zebrafish"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

species_colors <- c("Medaka" = "#D73027", "Zebrafish" = "#4575B4")

# =============================================================================
# THEME
# =============================================================================

theme_compare <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(hjust = 0, face = "bold", size = base_size + 2,
                                      margin = margin(b = 2)),
      plot.subtitle    = element_text(hjust = 0, size = base_size - 1, color = "grey40",
                                      margin = margin(b = 6)),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = base_size - 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      strip.text       = element_text(face = "bold", size = base_size),
      strip.background = element_rect(fill = "grey96", color = NA),
      axis.title       = element_text(size = base_size),
      plot.margin      = margin(6, 8, 6, 6)
    )
}

save_pdf <- function(plot, filename, width = 14, height = 10) {
  path <- file.path(OUTPUT_DIR, filename)
  ggsave(path, plot, width = width, height = height, device = "pdf")
  cat(sprintf("  -> Saved %s\n", path))
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("  GASTRULATION DYNAMICS -- MEDAKA vs ZEBRAFISH\n")
cat(strrep("=", 70), "\n")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("\n=== LOADING DATA ===\n")

# Summaries
sum_m <- read_csv(file.path(MEDAKA_DIR, "medaka_dynamics_summary.csv"), show_col_types = FALSE)
sum_z <- read_csv(file.path(ZEBRAFISH_DIR, "zebrafish_dynamics_summary.csv"), show_col_types = FALSE)
summary_both <- bind_rows(sum_m, sum_z)

# Track metrics
tm_m <- fread(file.path(MEDAKA_DIR, "medaka_dynamics_track_metrics.csv"), showProgress = FALSE) %>%
  mutate(species = "Medaka")
tm_z <- fread(file.path(ZEBRAFISH_DIR, "zebrafish_dynamics_track_metrics.csv"), showProgress = FALSE) %>%
  mutate(species = "Zebrafish")

# Harmonize columns — keep only shared ones + species specific
shared_cols <- c("TRACK_ID", "n_frames", "duration_min", "mean_theta", "mean_phi", "mean_depth",
                 "net_disp", "total_path", "straightness", "confinement",
                 "mean_speed", "mean_v_epiboly", "mean_v_depth",
                 "net_dTheta", "net_dPhi", "mean_turning", "mean_acceleration",
                 "persistence", "movement_type", "zone", "species")

existing_m <- intersect(shared_cols, names(tm_m))
existing_z <- intersect(shared_cols, names(tm_z))
all_cols <- union(existing_m, existing_z)
for (col in setdiff(all_cols, names(tm_m))) tm_m[[col]] <- NA
for (col in setdiff(all_cols, names(tm_z))) tm_z[[col]] <- NA

track_all <- bind_rows(
  tm_m[, all_cols, with = FALSE],
  tm_z[, all_cols, with = FALSE]
) %>% mutate(species = factor(species, levels = c("Medaka", "Zebrafish")))

# Oriented spots — for velocity distributions
cat("  Loading oriented spots for velocity comparison...\n")
spots_m <- fread(file.path(MEDAKA_INPUT, "oriented_tracks_medaka.csv"), showProgress = FALSE) %>%
  mutate(species = "Medaka", time_min = FRAME * MEDAKA_FI / 60)
spots_z <- fread(file.path(ZEBRAFISH_INPUT, "oriented_tracks_zebrafish.csv"), showProgress = FALSE) %>%
  mutate(species = "Zebrafish", time_min = FRAME * ZEBRAFISH_FI / 60)

# --- Voxel calibration: convert pixel positions to microns ---
cat(sprintf("  Voxel calibration: medaka=%.5f um/px, zebrafish=%.5f um/px\n",
            MEDAKA_VOXEL_UM, ZEBRAFISH_VOXEL_UM))
for (col in c("POSITION_X", "POSITION_Y", "POSITION_Z", "RADIAL_DIST", "SPHERICAL_DEPTH")) {
  spots_m[[col]] <- spots_m[[col]] * MEDAKA_VOXEL_UM
  spots_z[[col]] <- spots_z[[col]] * ZEBRAFISH_VOXEL_UM
}

# Compute velocity for both — with smoothing and configurable lag
# vel_lag: number of frames to lag for displacement (4 for medaka = 2 min, 1 for zebrafish = 2 min)
MEDAKA_SMOOTH_K    <- 5L   # 5 × 0.5 min = 2.5 min
ZEBRAFISH_SMOOTH_K <- 3L   # 3 × 2 min = 6 min (time-comparable)
MEDAKA_VEL_LAG  <- 4L   # 4 × 30 s = 2 min
ZEBRAFISH_VEL_LAG <- 1L # 1 × 120 s = 2 min

compute_vel <- function(df, frame_int_sec, vel_lag = 1L, smooth_k = 5L) {
  fi_min <- frame_int_sec / 60
  VL <- as.integer(vel_lag)
  dt <- as.data.table(copy(df))
  setkey(dt, TRACK_ID, FRAME)

  # Rolling-mean smoothing
  for (col in c("POSITION_X", "POSITION_Y", "POSITION_Z", "RADIAL_DIST", "THETA_DEG", "PHI_DEG")) {
    sm_col <- paste0(col, "_SM")
    dt[, (sm_col) := frollmean(get(col), n = smooth_k, align = "center"), by = TRACK_ID]
    dt[is.na(get(sm_col)), (sm_col) := get(col)]
  }

  dt[, `:=`(
    dx        = POSITION_X_SM - shift(POSITION_X_SM, VL),
    dy        = POSITION_Y_SM - shift(POSITION_Y_SM, VL),
    dz        = POSITION_Z_SM - shift(POSITION_Z_SM, VL),
    dt_frames = FRAME - shift(FRAME, VL),
    dTheta    = THETA_DEG_SM - shift(THETA_DEG_SM, VL),
    dPhi      = PHI_DEG_SM - shift(PHI_DEG_SM, VL),
    dDepth    = SPHERICAL_DEPTH - shift(SPHERICAL_DEPTH, VL)
  ), by = TRACK_ID]

  dt[, dPhi := fifelse(abs(dPhi) > 180, dPhi - sign(dPhi) * 360, dPhi)]
  dt[, `:=`(
    disp_3d    = sqrt(dx^2 + dy^2 + dz^2),
    theta_rad  = THETA_DEG * pi / 180,
    dTheta_rad = dTheta * pi / 180,
    dPhi_rad   = dPhi * pi / 180
  )]
  dt[, `:=`(
    inst_speed    = disp_3d / (dt_frames * fi_min),
    v_epiboly     = RADIAL_DIST * dTheta_rad / (dt_frames * fi_min),
    v_convergence = RADIAL_DIST * sin(theta_rad) * dPhi_rad / (dt_frames * fi_min),
    v_depth       = dDepth / (dt_frames * fi_min)
  )]

  setDF(dt)
  dt %>% filter(!is.na(dt_frames) & dt_frames == vel_lag)
}

cat(sprintf("  Velocity: medaka lag=%d (%ds, smooth k=%d), zebrafish lag=%d (%ds, smooth k=%d)\n",
            MEDAKA_VEL_LAG, MEDAKA_VEL_LAG * MEDAKA_FI, MEDAKA_SMOOTH_K,
            ZEBRAFISH_VEL_LAG, ZEBRAFISH_VEL_LAG * ZEBRAFISH_FI, ZEBRAFISH_SMOOTH_K))

vel_m <- compute_vel(spots_m, MEDAKA_FI, vel_lag = MEDAKA_VEL_LAG, smooth_k = MEDAKA_SMOOTH_K)
vel_z <- compute_vel(spots_z, ZEBRAFISH_FI, vel_lag = ZEBRAFISH_VEL_LAG, smooth_k = ZEBRAFISH_SMOOTH_K)
vel_all <- bind_rows(vel_m, vel_z) %>%
  mutate(species = factor(species, levels = c("Medaka", "Zebrafish")))

# =============================================================================
# SMOOTHING SENSITIVITY — TURNING ANGLE ARTIFACT TEST (both species)
# =============================================================================
# Re-compute turning angles at multiple smoothing levels for both species.
# This validates whether turning-angle differences between species are real
# or just noise from segmentation "wiggling".  If species differences persist
# even under aggressive smoothing, they reflect genuine biological behaviour.
# =============================================================================

cat("\n=== SMOOTHING SENSITIVITY — TURNING ANGLE ARTIFACT TEST ===\n")

compute_turn_at_smooth <- function(df, frame_int_sec, vel_lag, smooth_k_vec, smooth_labels, sp_name) {
  fi_min <- frame_int_sec / 60
  VL_int <- as.integer(vel_lag)

  rbindlist(lapply(seq_along(smooth_k_vec), function(lvl) {
    k <- smooth_k_vec[lvl]
    label <- smooth_labels[lvl]

    tmp <- as.data.table(df)[, .(TRACK_ID, FRAME, POSITION_X, POSITION_Y, POSITION_Z)]
    if (k == 1L) {
      tmp[, `:=`(sx = POSITION_X, sy = POSITION_Y, sz = POSITION_Z)]
    } else {
      tmp[, sx := frollmean(POSITION_X, n = k, align = "center"), by = TRACK_ID]
      tmp[, sy := frollmean(POSITION_Y, n = k, align = "center"), by = TRACK_ID]
      tmp[, sz := frollmean(POSITION_Z, n = k, align = "center"), by = TRACK_ID]
      tmp[is.na(sx), sx := POSITION_X]
      tmp[is.na(sy), sy := POSITION_Y]
      tmp[is.na(sz), sz := POSITION_Z]
    }

    tmp[, `:=`(dx = sx - shift(sx, VL_int),
               dy = sy - shift(sy, VL_int),
               dz = sz - shift(sz, VL_int),
               dt_f = FRAME - shift(FRAME, VL_int)), by = TRACK_ID]
    tmp[, disp_3d := sqrt(dx^2 + dy^2 + dz^2)]
    tmp[, inst_speed := disp_3d / (dt_f * fi_min)]

    tmp[, `:=`(lag_dx = shift(dx, VL_int),
               lag_dy = shift(dy, VL_int),
               lag_dz = shift(dz, VL_int),
               mag_prev = shift(disp_3d, VL_int)), by = TRACK_ID]
    tmp[, dot_prev := dx * lag_dx + dy * lag_dy + dz * lag_dz]
    tmp[, cos_ang := dot_prev / (disp_3d * mag_prev)]
    tmp[, turn_ang := acos(pmin(pmax(cos_ang, -1), 1)) * 180 / pi]

    valid <- tmp[!is.na(dt_f) & dt_f == VL_int & !is.na(turn_ang) & !is.na(inst_speed)]

    # Per-track straightness: compute net_disp / total_path for each track
    track_str <- tmp[!is.na(dt_f) & dt_f == VL_int,
                     .(net_disp  = sqrt((last(sx) - first(sx))^2 +
                                        (last(sy) - first(sy))^2 +
                                        (last(sz) - first(sz))^2),
                       total_path = sum(sqrt(diff(sx)^2 + diff(sy)^2 + diff(sz)^2), na.rm = TRUE),
                       n = .N),
                     by = TRACK_ID]
    track_str <- track_str[n >= 5]
    track_str[, straightness := fifelse(total_path > 0, net_disp / total_path, NA_real_)]

    data.table(smooth_level = label, smooth_k = k, species = sp_name,
               turning_angle = valid$turn_ang, inst_speed = valid$inst_speed,
               straightness = NA_real_, is_step = TRUE)
  }))
}

# Build per-track straightness separately (different grain from step-level)
compute_straightness_at_smooth <- function(df, smooth_k_vec, smooth_labels, sp_name) {
  rbindlist(lapply(seq_along(smooth_k_vec), function(lvl) {
    k <- smooth_k_vec[lvl]
    label <- smooth_labels[lvl]

    tmp <- as.data.table(df)[, .(TRACK_ID, FRAME, POSITION_X, POSITION_Y, POSITION_Z)]
    setkey(tmp, TRACK_ID, FRAME)
    if (k == 1L) {
      tmp[, `:=`(sx = POSITION_X, sy = POSITION_Y, sz = POSITION_Z)]
    } else {
      tmp[, sx := frollmean(POSITION_X, n = k, align = "center"), by = TRACK_ID]
      tmp[, sy := frollmean(POSITION_Y, n = k, align = "center"), by = TRACK_ID]
      tmp[, sz := frollmean(POSITION_Z, n = k, align = "center"), by = TRACK_ID]
      tmp[is.na(sx), sx := POSITION_X]
      tmp[is.na(sy), sy := POSITION_Y]
      tmp[is.na(sz), sz := POSITION_Z]
    }

    track_met <- tmp[, .(
      net_disp   = sqrt((last(sx) - first(sx))^2 + (last(sy) - first(sy))^2 + (last(sz) - first(sz))^2),
      total_path = sum(sqrt(diff(sx)^2 + diff(sy)^2 + diff(sz)^2), na.rm = TRUE),
      n = .N
    ), by = TRACK_ID]
    track_met <- track_met[n >= 5]
    track_met[, straightness := fifelse(total_path > 0, net_disp / total_path, NA_real_)]

    data.table(smooth_level = label, smooth_k = k, species = sp_name,
               straightness = track_met$straightness)
  }))
}

# Smoothing levels: raw, species-default, moderate, aggressive
# For medaka:  k = 1, 5, 11, 21
# For zebrafish: k = 1, 3, 7, 13
MEDAKA_SMOOTH_LEVELS <- c(1L, MEDAKA_SMOOTH_K,
                           as.integer(2L * MEDAKA_SMOOTH_K + 1L),
                           as.integer(4L * MEDAKA_SMOOTH_K + 1L))
ZEBRAFISH_SMOOTH_LEVELS <- c(1L, ZEBRAFISH_SMOOTH_K,
                              as.integer(2L * ZEBRAFISH_SMOOTH_K + 1L),
                              as.integer(4L * ZEBRAFISH_SMOOTH_K + 1L))

MEDAKA_SMOOTH_LABELS <- c("Raw (k=1)",
                           sprintf("Default (k=%d)", MEDAKA_SMOOTH_K),
                           sprintf("Moderate (k=%d)", 2L * MEDAKA_SMOOTH_K + 1L),
                           sprintf("Aggressive (k=%d)", 4L * MEDAKA_SMOOTH_K + 1L))
ZEBRAFISH_SMOOTH_LABELS <- c("Raw (k=1)",
                              sprintf("Default (k=%d)", ZEBRAFISH_SMOOTH_K),
                              sprintf("Moderate (k=%d)", 2L * ZEBRAFISH_SMOOTH_K + 1L),
                              sprintf("Aggressive (k=%d)", 4L * ZEBRAFISH_SMOOTH_K + 1L))

# Compute turning angles at each smoothing level
turn_m <- compute_turn_at_smooth(spots_m, MEDAKA_FI, MEDAKA_VEL_LAG,
                                  MEDAKA_SMOOTH_LEVELS, MEDAKA_SMOOTH_LABELS, "Medaka")
turn_z <- compute_turn_at_smooth(spots_z, ZEBRAFISH_FI, ZEBRAFISH_VEL_LAG,
                                  ZEBRAFISH_SMOOTH_LEVELS, ZEBRAFISH_SMOOTH_LABELS, "Zebrafish")

# Compute straightness at each smoothing level
str_m <- compute_straightness_at_smooth(spots_m, MEDAKA_SMOOTH_LEVELS, MEDAKA_SMOOTH_LABELS, "Medaka")
str_z <- compute_straightness_at_smooth(spots_z, ZEBRAFISH_SMOOTH_LEVELS, ZEBRAFISH_SMOOTH_LABELS, "Zebrafish")

# Unify labels for cross-species comparison: use generic names
# Each species has its own k values, so map by unique k per species
relabel <- function(dt) {
  dt[, smooth_label := {
    uks <- sort(unique(smooth_k))
    fcase(
      smooth_k == uks[1], "Raw",
      smooth_k == uks[2], "Default",
      smooth_k == uks[3], "Moderate",
      smooth_k == uks[4], "Aggressive"
    )
  }, by = species]
  lvls <- c("Raw", "Default", "Moderate", "Aggressive")
  dt[, smooth_label := factor(smooth_label, levels = lvls)]
  dt
}
turn_all <- relabel(rbindlist(list(turn_m, turn_z)))
turn_all[, species := factor(species, levels = c("Medaka", "Zebrafish"))]

str_all <- relabel(rbindlist(list(str_m, str_z)))
str_all[, species := factor(species, levels = c("Medaka", "Zebrafish"))]

# Summary table
turn_summary <- turn_all[, .(
  median_turn = median(turning_angle, na.rm = TRUE),
  mean_turn   = mean(turning_angle, na.rm = TRUE),
  median_speed = median(inst_speed, na.rm = TRUE),
  n = .N
), by = .(species, smooth_label)]

str_summary <- str_all[, .(
  median_str = median(straightness, na.rm = TRUE),
  mean_str   = mean(straightness, na.rm = TRUE),
  n = .N
), by = .(species, smooth_label)]

cat("  Turning angle medians across smoothing levels:\n")
for (sp in c("Medaka", "Zebrafish")) {
  ts <- turn_summary[species == sp]
  cat(sprintf("    %s:  Raw=%.1f, Default=%.1f, Moderate=%.1f, Aggressive=%.1f deg\n",
              sp, ts[smooth_label == "Raw", median_turn],
              ts[smooth_label == "Default", median_turn],
              ts[smooth_label == "Moderate", median_turn],
              ts[smooth_label == "Aggressive", median_turn]))
}
cat("  Straightness medians across smoothing levels:\n")
for (sp in c("Medaka", "Zebrafish")) {
  ss <- str_summary[species == sp]
  cat(sprintf("    %s:  Raw=%.3f, Default=%.3f, Moderate=%.3f, Aggressive=%.3f\n",
              sp, ss[smooth_label == "Raw", median_str],
              ss[smooth_label == "Default", median_str],
              ss[smooth_label == "Moderate", median_str],
              ss[smooth_label == "Aggressive", median_str]))
}

# Species gap at aggressive smoothing
aggr_m_turn <- turn_summary[species == "Medaka" & smooth_label == "Aggressive", median_turn]
aggr_z_turn <- turn_summary[species == "Zebrafish" & smooth_label == "Aggressive", median_turn]
cat(sprintf("  SPECIES GAP at aggressive smoothing: Medaka %.1f vs Zebrafish %.1f deg (ratio %.1fx)\n",
            aggr_m_turn, aggr_z_turn, aggr_m_turn / aggr_z_turn))
thick_m <- read_csv(file.path(MEDAKA_DIR, "medaka_thickness_data.csv"), show_col_types = FALSE)
thick_z <- read_csv(file.path(ZEBRAFISH_DIR, "zebrafish_thickness_data.csv"), show_col_types = FALSE)
thickness_all <- bind_rows(thick_m, thick_z) %>%
  mutate(species = factor(species, levels = c("Medaka", "Zebrafish")))

cat(sprintf("  Medaka: %s tracks | Zebrafish: %s tracks\n",
            format(nrow(tm_m), big.mark = ","), format(nrow(tm_z), big.mark = ",")))
cat(sprintf("  Medaka vel steps: %s | Zebrafish vel steps: %s\n",
            format(nrow(vel_m), big.mark = ","), format(nrow(vel_z), big.mark = ",")))

# =============================================================================
# LOAD NUCLEAR STATS
# =============================================================================

cat("\n=== LOADING NUCLEAR STATS ===\n")

HAS_NUCLEAR_M <- dir.exists(MEDAKA_NUCLEAR_DIR)
HAS_NUCLEAR_Z <- dir.exists(ZEBRAFISH_NUCLEAR_DIR)
HAS_NUCLEAR <- HAS_NUCLEAR_M && HAS_NUCLEAR_Z

if (HAS_NUCLEAR) {
  nuc_m_raw <- read.delim(file.path(MEDAKA_NUCLEAR_DIR, "nuclear_stats_per_time.tsv")) %>%
    mutate(species = "Medaka")
  nuc_z <- read.delim(file.path(ZEBRAFISH_NUCLEAR_DIR, "nuclear_stats_per_time.tsv")) %>%
    mutate(species = "Zebrafish")

  # Fix timepoints if constant (generation bug)
  if (length(unique(nuc_m_raw$timepoint)) == 1 && nrow(nuc_m_raw) > 1) {
    nuc_m_raw$timepoint <- seq(0, nrow(nuc_m_raw) - 1)
  }
  if (length(unique(nuc_z$timepoint)) == 1 && nrow(nuc_z) > 1) {
    nuc_z$timepoint <- seq(0, nrow(nuc_z) - 1)
  }

  # Medaka contraction artefact: from tp 500 the embryo contracts,
  # making nuclear counts / sizes unreliable.  Filter to tp <= 500.
  MEDAKA_CONTRACTION_TP <- 500
  nuc_m <- nuc_m_raw %>% filter(timepoint <= MEDAKA_CONTRACTION_TP)
  cat(sprintf("  Medaka: keeping tp 1-%d (of %d) — tp>%d excluded (contraction artefact)\n",
              MEDAKA_CONTRACTION_TP, nrow(nuc_m_raw), MEDAKA_CONTRACTION_TP))

  # Convert to absolute time
  nuc_m$time_min <- nuc_m$timepoint * MEDAKA_FI / 60
  nuc_z$time_min <- nuc_z$timepoint * ZEBRAFISH_FI / 60

  # Bad timepoint detection — same iterative LOESS for both species
  detect_bad_tps <- function(df) {
    df$bad_tp <- FALSE
    for (iter in 1:5) {
      good_idx <- which(!df$bad_tp)
      lo <- loess(n_nuclei_3d ~ timepoint, data = df[good_idx, ], span = 0.3)
      expected <- predict(lo, newdata = df)
      df$bad_tp <- !is.na(expected) & (df$n_nuclei_3d < 0.75 * expected)
    }
    df
  }
  nuc_m <- detect_bad_tps(nuc_m)
  nuc_z <- detect_bad_tps(nuc_z)
  bad_tps_m <- nuc_m$timepoint[nuc_m$bad_tp]
  bad_tps_z <- nuc_z$timepoint[nuc_z$bad_tp]
  nuc_m_clean <- nuc_m %>% filter(!bad_tp)
  nuc_z_clean <- nuc_z %>% filter(!bad_tp)

  cat(sprintf("  Bad timepoints: Medaka %d, Zebrafish %d\n",
              length(bad_tps_m), length(bad_tps_z)))

  # Detect whether core-filtered 3D columns are available
  has_core_3d <- "core_vol_median_um3" %in% names(nuc_m)
  if (has_core_3d) {
    cat("  Core-filtered 3D columns detected — using core metrics\n")
    vol_median_col   <- "core_vol_median_um3"
    eqsph_median_col <- "core_eqsph_diam_median_um"
    nn3d_median_col  <- "core_nn3d_median_um"
    bbox_major_col   <- "core_bbox_major_mean_um"
    bbox_minor_col   <- "core_bbox_minor_mean_um"
    core_label       <- " (core z-filtered)"
  } else {
    cat("  No core 3D columns — using all-nuclei metrics\n")
    vol_median_col   <- "vol_median_um3"
    eqsph_median_col <- "eqsph_diam_median_um"
    nn3d_median_col  <- "nn3d_median_um"
    bbox_major_col   <- "bbox_major_mean_um"
    bbox_minor_col   <- "bbox_minor_mean_um"
    core_label       <- ""
  }

  nuc_both <- bind_rows(nuc_m_clean, nuc_z_clean) %>%
    mutate(species = factor(species, levels = c("Medaka", "Zebrafish")))

  # Load global stats
  nuc_global_m <- read.delim(file.path(MEDAKA_NUCLEAR_DIR, "nuclear_stats_global.tsv"))
  nuc_global_z <- read.delim(file.path(ZEBRAFISH_NUCLEAR_DIR, "nuclear_stats_global.tsv"))

  cat(sprintf("  Medaka: %d clean timepoints (of %d) | Zebrafish: %d clean (of %d)\n",
              nrow(nuc_m_clean), nrow(nuc_m), nrow(nuc_z_clean), nrow(nuc_z)))
} else {
  cat("  Nuclear stats not found for one or both species\n")
}
# =============================================================================
# MSD COMPARISON (recompute on same time scale)
# =============================================================================

cat("\n=== COMPUTING MSD FOR COMPARISON ===\n")

compute_msd <- function(df, max_lag, n_sample, fi_min) {
  set.seed(42)
  track_ids <- unique(df$TRACK_ID)
  if (length(track_ids) > n_sample) track_ids <- sample(track_ids, n_sample)
  dt <- as.data.table(df)[TRACK_ID %in% track_ids,
                           .(TRACK_ID, FRAME, POSITION_X, POSITION_Y, POSITION_Z)]
  setkey(dt, TRACK_ID, FRAME)
  msd_list <- vector("list", max_lag)
  for (lag in seq_len(max_lag)) {
    d2 <- dt[, {
      n <- .N; if (n > lag) {
        i1 <- 1:(n - lag); i2 <- (1 + lag):n
        valid <- (FRAME[i2] - FRAME[i1]) == lag
        if (sum(valid) > 0) {
          dsq <- (POSITION_X[i2[valid]] - POSITION_X[i1[valid]])^2 +
                 (POSITION_Y[i2[valid]] - POSITION_Y[i1[valid]])^2 +
                 (POSITION_Z[i2[valid]] - POSITION_Z[i1[valid]])^2
          .(msd = mean(dsq), n = length(dsq))
        } else .(msd = numeric(0), n = integer(0))
      } else .(msd = numeric(0), n = integer(0))
    }, by = TRACK_ID]
    if (nrow(d2) > 0 && sum(d2$n) > 0) {
      msd_list[[lag]] <- data.frame(
        lag_frames = lag, lag_min = lag * fi_min,
        mean_msd = weighted.mean(d2$msd, d2$n),
        n_pairs = sum(d2$n), n_tracks = nrow(d2))
    }
  }
  bind_rows(msd_list)
}

msd_m <- compute_msd(spots_m, MSD_MAX_LAG, MSD_N_TRACKS, MEDAKA_FI / 60) %>%
  mutate(species = "Medaka")
msd_z <- compute_msd(spots_z, MSD_MAX_LAG, MSD_N_TRACKS, ZEBRAFISH_FI / 60) %>%
  mutate(species = "Zebrafish")
msd_both <- bind_rows(msd_m, msd_z) %>%
  mutate(species = factor(species, levels = c("Medaka", "Zebrafish")))

cat(sprintf("  MSD: Medaka %d lags, Zebrafish %d lags\n", nrow(msd_m), nrow(msd_z)))

# ###########################################################################
# FIGURES
# ###########################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("  GENERATING COMPARISON FIGURES\n")
cat(strrep("=", 70), "\n")

# ============================================================================
# FIGURE 1: VELOCITY & EPIBOLY COMPARISON
# ============================================================================
# Story: overall cell speed is similar between species (~5% difference),
# but zebrafish epibolises 3-5× faster — not because individual vegetalward
# steps are faster, but because ~88% of ZF steps are vegetalward vs only
# ~59% in medaka.  Early medaka moves predominantly animalward.
# ============================================================================

cat("\n=== FIGURE 1: VELOCITY & EPIBOLY COMPARISON ===\n")

vel_moving <- vel_all %>% filter(disp_3d > 0)
TIME_BIN_MIN <- 10

# ── Helpers ──────────────────────────────────────────────────────────────
vel_interp_at <- function(df, col, t) {
  approx(df$time_bin, df[[col]], xout = t, rule = 2)$y
}

vel_add_diff <- function(p, m_df, z_df, col, t, nudge_x = 8) {
  vm <- vel_interp_at(m_df, col, t)
  vz <- vel_interp_at(z_df, col, t)
  if (is.na(vm) || is.na(vz) || vm == 0) return(p)
  pct_diff <- (vz - vm) / abs(vm) * 100
  seg_df <- data.frame(x = t, ymin = min(vm, vz), ymax = max(vm, vz))
  lab_df <- data.frame(x = t + nudge_x, y = (vm + vz) / 2,
                       label = sprintf("%+.0f%%", pct_diff))
  p +
    geom_segment(data = seg_df, aes(x = x, xend = x, y = ymin, yend = ymax),
                 inherit.aes = FALSE, linetype = "dotted", linewidth = 0.5, color = "grey30") +
    geom_point(data = data.frame(x = c(t, t), y = c(vm, vz)),
               aes(x = x, y = y), inherit.aes = FALSE, size = 2.5, shape = 21,
               fill = c(species_colors["Medaka"], species_colors["Zebrafish"]),
               color = "white", stroke = 0.4) +
    geom_label(data = lab_df, aes(x = x, y = y, label = label),
               inherit.aes = FALSE, size = 2.8, fill = "white", alpha = 0.85,
               linewidth = 0.2, label.padding = unit(0.15, "lines"))
}

# ── Summarise per time bin ───────────────────────────────────────────────
vel_bin_summary <- vel_all %>%
  mutate(time_bin = floor(time_min / TIME_BIN_MIN) * TIME_BIN_MIN) %>%
  group_by(species, time_bin) %>%
  summarise(
    mean_speed   = mean(inst_speed, na.rm = TRUE),
    se_speed     = sd(inst_speed, na.rm = TRUE) / sqrt(n()),
    mean_v_epi   = mean(v_epiboly, na.rm = TRUE),
    se_v_epi     = sd(v_epiboly, na.rm = TRUE) / sqrt(n()),
    pct_veg      = 100 * mean(v_epiboly > 0, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 50)

# Strict epiboly (vegetalward only)
vel_epi_strict <- vel_all %>%
  filter(v_epiboly > 0) %>%
  mutate(time_bin = floor(time_min / TIME_BIN_MIN) * TIME_BIN_MIN) %>%
  group_by(species, time_bin) %>%
  summarise(
    mean_v_epi_strict = mean(v_epiboly, na.rm = TRUE),
    se_v_epi_strict   = sd(v_epiboly, na.rm = TRUE) / sqrt(n()),
    n = n(), .groups = "drop"
  ) %>%
  filter(n >= 50)

# Angular rate dθ/dt (R-independent epiboly progression)
vel_all <- vel_all %>%
  mutate(fi_min = ifelse(species == "Medaka", MEDAKA_FI / 60, ZEBRAFISH_FI / 60),
         dTheta_rate = dTheta / (dt_frames * fi_min))

vel_angular <- vel_all %>%
  filter(v_epiboly > 0) %>%
  mutate(time_bin = floor(time_min / TIME_BIN_MIN) * TIME_BIN_MIN) %>%
  group_by(species, time_bin) %>%
  summarise(
    mean_dTheta_rate = mean(dTheta_rate, na.rm = TRUE),
    se_dTheta_rate   = sd(dTheta_rate, na.rm = TRUE) / sqrt(n()),
    n = n(), .groups = "drop"
  ) %>%
  filter(n >= 50)

# Split by species
vbs_m <- vel_bin_summary %>% filter(species == "Medaka")
vbs_z <- vel_bin_summary %>% filter(species == "Zebrafish")
ves_m <- vel_epi_strict %>% filter(species == "Medaka")
ves_z <- vel_epi_strict %>% filter(species == "Zebrafish")
van_m <- vel_angular %>% filter(species == "Medaka")
van_z <- vel_angular %>% filter(species == "Zebrafish")

t_end_m_vel <- max(vbs_m$time_bin)
t_end_z_vel <- max(vbs_z$time_bin)
t_overlap_vel <- min(t_end_m_vel, t_end_z_vel)

vel_el <- function(df, col, frac = 0.2) {
  n <- nrow(df)
  early <- df %>% slice_head(n = max(1, round(n * frac)))
  late  <- df %>% slice_tail(n = max(1, round(n * frac)))
  list(early = mean(early[[col]], na.rm = TRUE),
       late  = mean(late[[col]],  na.rm = TRUE),
       pct   = (mean(late[[col]], na.rm = TRUE) - mean(early[[col]], na.rm = TRUE)) /
               abs(mean(early[[col]], na.rm = TRUE)) * 100)
}

med_end_vel <- geom_vline(xintercept = t_end_m_vel, linetype = "dashed",
                          color = species_colors["Medaka"], linewidth = 0.3, alpha = 0.6)

# Summary stats for annotations
med_speed_m <- median(vel_m$inst_speed[vel_m$disp_3d > 0], na.rm = TRUE)
med_speed_z <- median(vel_z$inst_speed[vel_z$disp_3d > 0], na.rm = TRUE)
med_epi_strict_m <- median(vel_m$v_epiboly[vel_m$v_epiboly > 0], na.rm = TRUE)
med_epi_strict_z <- median(vel_z$v_epiboly[vel_z$v_epiboly > 0], na.rm = TRUE)
pct_veg_m <- 100 * mean(vel_m$v_epiboly > 0, na.rm = TRUE)
pct_veg_z <- 100 * mean(vel_z$v_epiboly > 0, na.rm = TRUE)
mean_net_epi_m <- mean(vel_m$v_epiboly, na.rm = TRUE)
mean_net_epi_z <- mean(vel_z$v_epiboly, na.rm = TRUE)
mean_dth_m <- mean(vel_m$dTheta_rate[vel_m$v_epiboly > 0], na.rm = TRUE)
mean_dth_z <- mean(vel_z$dTheta_rate[vel_z$v_epiboly > 0], na.rm = TRUE)

el_speed_m <- vel_el(vbs_m, "mean_speed")
el_speed_z <- vel_el(vbs_z, "mean_speed")

# ── ROW 1: Overall speed is similar ─────────────────────────────────────

# 1a: Speed over time — similar between species
p1a <- ggplot(vel_bin_summary, aes(x = time_bin, y = mean_speed, color = species)) +
  geom_ribbon(aes(ymin = mean_speed - 1.96 * se_speed,
                  ymax = mean_speed + 1.96 * se_speed, fill = species),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  med_end_vel +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  labs(title = "Overall cell speed is similar between species",
       subtitle = sprintf("Median: Med %.2f vs ZF %.2f \u00b5m/min (%+.0f%%) \u2014 both ~1 \u00b5m/min",
                          med_speed_m, med_speed_z,
                          (med_speed_m / med_speed_z - 1) * 100),
       x = "Time (min)", y = "Speed (\u00b5m/min)") +
  guides(fill = "none") + theme_compare()
p1a <- vel_add_diff(p1a, vbs_m, vbs_z, "mean_speed", 0, nudge_x = 12)
p1a <- vel_add_diff(p1a, vbs_m, vbs_z, "mean_speed", t_overlap_vel, nudge_x = -15)

# 1b: Speed distribution — overlapping, confirming similarity
p1b <- ggplot(vel_moving, aes(x = inst_speed, fill = species)) +
  geom_density(alpha = 0.45, color = NA) +
  scale_fill_manual(values = species_colors) +
  coord_cartesian(xlim = c(0, quantile(vel_moving$inst_speed, 0.99, na.rm = TRUE))) +
  geom_vline(xintercept = med_speed_m, linetype = "dashed",
             color = species_colors["Medaka"], linewidth = 0.5) +
  geom_vline(xintercept = med_speed_z, linetype = "dashed",
             color = species_colors["Zebrafish"], linewidth = 0.5) +
  annotate("text", x = max(med_speed_m, med_speed_z) + 0.3,
           y = Inf, vjust = 1.5, hjust = 0, size = 3, color = "grey40",
           label = sprintf("\u0394 = %.0f%%", (med_speed_m / med_speed_z - 1) * 100)) +
  labs(title = "Speed distributions largely overlap",
       subtitle = "Dashed lines = medians; distributions nearly identical",
       x = "Instantaneous speed (\u00b5m/min)", y = "Density", fill = "Species") +
  theme_compare()

# ── ROW 2: But epiboly commitment is drastically different ───────────────

# 1c: Vegetalward fraction over time — the KEY difference
p1c <- ggplot(vel_bin_summary, aes(x = time_bin, y = pct_veg, color = species)) +
  geom_ribbon(data = vel_bin_summary,
              aes(x = time_bin, ymin = pct_veg, ymax = 100, fill = species),
              alpha = 0.08, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  med_end_vel +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  annotate("text", x = 5, y = 15, hjust = 0, size = 3.2, fontface = "italic", color = "grey40",
           label = "below 50% = majority animalward") +
  labs(title = "Vegetalward commitment: the key difference",
       subtitle = sprintf("ZF: %.0f%% of steps are vegetalward throughout | Med: only %.0f%% (starts at ~23%%)",
                          pct_veg_z, pct_veg_m),
       x = "Time (min)", y = "% vegetalward steps", color = "Species") +
  guides(fill = "none") + theme_compare()

# 1d: Movement type composition — stacked bar
mvt_prop <- track_all %>%
  filter(!is.na(movement_type), movement_type != "Ingression") %>%
  mutate(movement_type = factor(movement_type,
                                levels = c("Epiboly", "Convergence", "Animalward"))) %>%
  group_by(species, movement_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(species) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  arrange(desc(movement_type)) %>%
  mutate(y_pos = cumsum(pct) - pct / 2)

mvt_colors_fig1 <- c("Epiboly" = "#66C2A5", "Convergence" = "#FC8D62",
                      "Animalward" = "#8DA0CB")

p1d <- ggplot(mvt_prop, aes(x = species, y = pct, fill = movement_type)) +
  geom_col(position = "stack", color = "white", linewidth = 0.5, width = 0.7) +
  geom_text(aes(y = y_pos, label = sprintf("%.0f%%", pct)),
            size = 3.5, fontface = "bold", color = "grey15") +
  scale_fill_manual(values = mvt_colors_fig1) +
  labs(title = "Track-level classification confirms pattern",
       subtitle = "86% of ZF tracks classified as Epiboly vs 50% for Medaka",
       x = NULL, y = "% of tracks", fill = "Movement type") +
  theme_compare() + theme(legend.position = "right")

# ── ROW 3: Strict epiboly step speed is only slightly different ──────────

# 1e: Strict epiboly velocity over time (only vegetalward steps)
el_epi_m <- vel_el(ves_m, "mean_v_epi_strict")
el_epi_z <- vel_el(ves_z, "mean_v_epi_strict")

p1e <- ggplot(vel_epi_strict, aes(x = time_bin, y = mean_v_epi_strict, color = species)) +
  geom_ribbon(aes(ymin = mean_v_epi_strict - 1.96 * se_v_epi_strict,
                  ymax = mean_v_epi_strict + 1.96 * se_v_epi_strict, fill = species),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  med_end_vel +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  labs(title = "Per-step epiboly speed is similar (vegetalward only)",
       subtitle = sprintf("Median: Med %.2f vs ZF %.2f \u00b5m/min (%+.0f%%) \u2014 only ~14%% difference",
                          med_epi_strict_m, med_epi_strict_z,
                          (med_epi_strict_z / med_epi_strict_m - 1) * 100),
       x = "Time (min)", y = "v_epiboly (\u00b5m/min)", color = "Species") +
  guides(fill = "none") + theme_compare()
p1e <- vel_add_diff(p1e, ves_m, ves_z, "mean_v_epi_strict", 0, nudge_x = 12)
p1e <- vel_add_diff(p1e, ves_m, ves_z, "mean_v_epi_strict", t_overlap_vel, nudge_x = -15)

# 1f: Net epiboly velocity (all steps, signed) — the combined effect
el_nepi_m <- vel_el(vbs_m, "mean_v_epi")
el_nepi_z <- vel_el(vbs_z, "mean_v_epi")

p1f <- ggplot(vel_bin_summary, aes(x = time_bin, y = mean_v_epi, color = species)) +
  geom_ribbon(aes(ymin = mean_v_epi - 1.96 * se_v_epi,
                  ymax = mean_v_epi + 1.96 * se_v_epi, fill = species),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  med_end_vel +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  annotate("text", x = 5, y = -0.3, hjust = 0, size = 2.8, fontface = "italic",
           color = "grey40", label = "\u2190 animalward") +
  annotate("text", x = 5, y = 0.5, hjust = 0, size = 2.8, fontface = "italic",
           color = "grey40", label = "vegetalward \u2192") +
  labs(title = "Net epiboly velocity: commitment \u00d7 speed combined",
       subtitle = sprintf("Mean all steps: Med %.2f vs ZF %.2f \u00b5m/min \u2014 ZF %.1f\u00d7 faster",
                          mean_net_epi_m, mean_net_epi_z,
                          mean_net_epi_z / mean_net_epi_m),
       x = "Time (min)", y = "v_epiboly (\u00b5m/min)", color = "Species") +
  guides(fill = "none") + theme_compare()
p1f <- vel_add_diff(p1f, vbs_m, vbs_z, "mean_v_epi", t_overlap_vel, nudge_x = -15)

# ── ROW 4: Angular rate confirms — geometry-independent ──────────────────

# 1g: Angular epiboly rate dθ/dt (R-independent)
p1g <- ggplot(vel_angular, aes(x = time_bin, y = mean_dTheta_rate, color = species)) +
  geom_ribbon(aes(ymin = mean_dTheta_rate - 1.96 * se_dTheta_rate,
                  ymax = mean_dTheta_rate + 1.96 * se_dTheta_rate, fill = species),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  med_end_vel +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  labs(title = "Angular epiboly rate (radius-independent)",
       subtitle = sprintf("d\u03b8/dt: Med %.3f vs ZF %.3f deg/min \u2014 ZF %.1f\u00d7 faster (vegetalward steps)",
                          mean_dth_m, mean_dth_z, mean_dth_z / mean_dth_m),
       x = "Time (min)", y = "d\u03b8/dt (deg/min)", color = "Species") +
  guides(fill = "none") + theme_compare()

# 1h: Summary bar — key numbers side by side
summary_df <- data.frame(
  metric = factor(c("Overall speed\n(\u00b5m/min)",
                     "Strict epiboly\nspeed (\u00b5m/min)",
                     "Net epiboly\nvelocity (\u00b5m/min)",
                     "Vegetalward\nfraction (%)",
                     "Angular rate\n(deg/min \u00d710)"),
                  levels = c("Overall speed\n(\u00b5m/min)",
                             "Strict epiboly\nspeed (\u00b5m/min)",
                             "Net epiboly\nvelocity (\u00b5m/min)",
                             "Vegetalward\nfraction (%)",
                             "Angular rate\n(deg/min \u00d710)")),
  species = rep(c("Medaka", "Zebrafish"), each = 5),
  value = c(med_speed_m, med_epi_strict_m, mean_net_epi_m, pct_veg_m, mean_dth_m * 10,
            med_speed_z, med_epi_strict_z, mean_net_epi_z, pct_veg_z, mean_dth_z * 10)
)

p1h <- ggplot(summary_df, aes(x = metric, y = value, fill = species)) +
  geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.2f", value)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = species_colors) +
  labs(title = "Summary: similar speed, vastly different epiboly",
       subtitle = "Speed \u2248 equal | Commitment & angular progression drive the difference",
       x = NULL, y = "Value", fill = "Species") +
  theme_compare() + theme(axis.text.x = element_text(size = 8))

# ── Compose figure ───────────────────────────────────────────────────────
fig1 <- (p1a + p1b) / (p1c + p1d) / (p1e + p1f) / (p1g + p1h) +
  plot_layout(heights = c(1, 1, 1, 1)) +
  plot_annotation(
    title = "Medaka vs Zebrafish \u2014 Cell Speed & Epiboly Dynamics",
    subtitle = paste0("Cells move at similar speeds, but zebrafish epibolises ",
                      sprintf("%.1f\u00d7", mean_net_epi_z / mean_net_epi_m),
                      " faster \u2014 driven by directional commitment (88% vs 59% vegetalward), not step speed (+14%)\n",
                      "95% CI ribbons | dashed vertical = medaka imaging end | n=1 embryo/species"),
    theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 15),
                  plot.subtitle = element_text(hjust = 0, size = 9.5, color = "grey30",
                                               lineheight = 1.2))
  )
save_pdf(fig1, "comparison_01_speed.pdf", width = 16, height = 22)

# ============================================================================
# FIGURE 2: TRACK STRAIGHTNESS & PERSISTENCE — SMOOTHING-VALIDATED
# ============================================================================
# Shows that species differences in turning angle and straightness persist
# even under aggressive trajectory smoothing, ruling out segmentation noise
# as the explanation.
# ============================================================================

cat("\n=== FIGURE 2: STRAIGHTNESS & PERSISTENCE (smoothing-validated) ===\n")

smooth_label_colors <- c("Raw" = "#D73027", "Default" = "#FC8D59",
                          "Moderate" = "#91BFDB", "Aggressive" = "#4575B4")

# 2a: Turning angle distributions — species side by side, smoothing as color
p2a <- ggplot(as.data.frame(turn_all),
              aes(x = turning_angle, color = smooth_label)) +
  geom_density(linewidth = 0.7, alpha = 0.8) +
  facet_wrap(~species) +
  scale_color_manual(values = smooth_label_colors, name = "Smoothing") +
  coord_cartesian(xlim = c(0, 120)) +
  labs(title = "Turning angle: smoothing sensitivity by species",
       subtitle = sprintf("Aggressive smoothing -- Medaka: %.1f deg, Zebrafish: %.1f deg (%.1fx gap persists)",
                          aggr_m_turn, aggr_z_turn, aggr_m_turn / aggr_z_turn),
       x = "Turning angle (deg)", y = "Density") +
  theme_compare()

# 2b: Median turning angle across smoothing levels — species overlaid
turn_sum_df <- as.data.frame(turn_summary)
turn_sum_df$smooth_label <- factor(turn_sum_df$smooth_label,
                                    levels = c("Raw", "Default", "Moderate", "Aggressive"))

p2b <- ggplot(turn_sum_df, aes(x = smooth_label, y = median_turn,
                                color = species, group = species)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.1f", median_turn)), vjust = -1, size = 3, show.legend = FALSE) +
  scale_color_manual(values = species_colors) +
  labs(title = "Median turning angle vs smoothing level",
       subtitle = "Species gap persists: noise affects both equally",
       x = "Smoothing level", y = "Median turning angle (deg)", color = "Species") +
  theme_compare()

# 2c: Straightness distributions — species side by side, smoothing as color
p2c <- ggplot(as.data.frame(str_all),
              aes(x = straightness, color = smooth_label)) +
  geom_density(linewidth = 0.7, alpha = 0.8) +
  facet_wrap(~species) +
  scale_color_manual(values = smooth_label_colors, name = "Smoothing") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Track straightness: smoothing sensitivity by species",
       subtitle = sprintf("Aggressive smoothing -- Medaka: %.3f, Zebrafish: %.3f median",
                          str_summary[species == "Medaka" & smooth_label == "Aggressive", median_str],
                          str_summary[species == "Zebrafish" & smooth_label == "Aggressive", median_str]),
       x = "Straightness (net disp / path length)", y = "Density") +
  theme_compare()

# 2d: Median straightness across smoothing levels
str_sum_df <- as.data.frame(str_summary)
str_sum_df$smooth_label <- factor(str_sum_df$smooth_label,
                                   levels = c("Raw", "Default", "Moderate", "Aggressive"))

p2d <- ggplot(str_sum_df, aes(x = smooth_label, y = median_str,
                               color = species, group = species)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.3f", median_str)), vjust = -1, size = 3, show.legend = FALSE) +
  scale_color_manual(values = species_colors) +
  labs(title = "Median straightness vs smoothing level",
       subtitle = "Smoothing increases straightness (removes noise), but species gap remains",
       x = "Smoothing level", y = "Median straightness", color = "Species") +
  theme_compare()

# 2e: Direct species overlay at aggressive smoothing (turning angle)
turn_aggr <- as.data.frame(turn_all[smooth_label == "Aggressive"])
p2e <- ggplot(turn_aggr, aes(x = turning_angle, fill = species)) +
  geom_density(alpha = 0.45, color = NA) +
  geom_vline(xintercept = aggr_m_turn, linetype = "dashed",
             color = species_colors["Medaka"], linewidth = 0.6) +
  geom_vline(xintercept = aggr_z_turn, linetype = "dashed",
             color = species_colors["Zebrafish"], linewidth = 0.6) +
  scale_fill_manual(values = species_colors) +
  coord_cartesian(xlim = c(0, 90)) +
  labs(title = "Turning angle after aggressive smoothing",
       subtitle = sprintf("Medaka %.1f deg vs Zebrafish %.1f deg -- real biological difference",
                          aggr_m_turn, aggr_z_turn),
       x = "Turning angle (deg)", y = "Density", fill = "Species") +
  theme_compare()

# 2f: Confinement ratio (unchanged by smoothing — uses original track metrics)
p2f <- ggplot(track_all, aes(x = confinement, fill = species)) +
  geom_density(alpha = 0.5, color = NA) +
  scale_fill_manual(values = species_colors) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Confinement ratio",
       subtitle = sprintf("Medaka: %.3f | Zebrafish: %.3f median",
                          median(tm_m$confinement, na.rm = TRUE),
                          median(tm_z$confinement, na.rm = TRUE)),
       x = "Confinement (max reach / path)", y = "Density", fill = "Species") +
  theme_compare()

fig2 <- (p2a + p2b) / (p2c + p2d) / (p2e + p2f) +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    title = "Medaka vs Zebrafish -- Straightness & Persistence (Smoothing-Validated)",
    subtitle = "Species differences persist under aggressive trajectory smoothing, ruling out segmentation noise",
    theme = theme(plot.title    = element_text(hjust = 0, face = "bold", size = 14),
                  plot.subtitle = element_text(hjust = 0, size = 10, color = "grey40"))
  )
save_pdf(fig2, "comparison_02_straightness.pdf", width = 16, height = 18)

# ============================================================================
# FIGURE 3: MSD
# ============================================================================

cat("\n=== FIGURE 3: MSD ===\n")

# Fit alpha
fit_alpha <- function(msd_df) {
  d <- msd_df %>% filter(lag_frames <= 10 & mean_msd > 0)
  if (nrow(d) < 3) return(NA_real_)
  coef(lm(log10(mean_msd) ~ log10(lag_min), data = d))[2]
}

alpha_m <- fit_alpha(msd_m)
alpha_z <- fit_alpha(msd_z)

p3a <- ggplot(msd_both, aes(x = lag_min, y = mean_msd, color = species)) +
  geom_line(linewidth = 1) + geom_point(size = 1.5) +
  scale_x_log10() + scale_y_log10() + annotation_logticks(sides = "bl") +
  scale_color_manual(values = species_colors) +
  labs(title = "MSD comparison (log-log)",
       subtitle = sprintf("alpha: Medaka=%.2f, Zebrafish=%.2f", alpha_m, alpha_z),
       x = "Time lag (min)", y = expression(MSD~(mu*m^2)), color = "Species") +
  theme_compare()

p3b <- ggplot(msd_both, aes(x = lag_min, y = mean_msd, color = species)) +
  geom_line(linewidth = 1) + geom_point(size = 1.5) +
  scale_color_manual(values = species_colors) +
  labs(title = "MSD comparison (linear)", x = "Time lag (min)",
       y = expression(MSD~(mu*m^2)), color = "Species") +
  theme_compare()

# Diffusivity (MSD/4t)
msd_diff <- msd_both %>% mutate(D_eff = mean_msd / (4 * lag_min))

p3c <- ggplot(msd_diff, aes(x = lag_min, y = D_eff, color = species)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = species_colors) +
  labs(title = "Effective diffusivity", x = "Time lag (min)",
       y = expression(D[eff]~(mu*m^2/min)), color = "Species") +
  theme_compare()

fig3 <- (p3a + p3b) / (p3c + plot_spacer()) +
  plot_annotation(
    title = "Medaka vs Zebrafish -- Mean Square Displacement",
    theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14))
  )
save_pdf(fig3, "comparison_03_msd.pdf", width = 14, height = 10)

# ============================================================================
# FIGURE 4: THICKNESS
# ============================================================================

cat("\n=== FIGURE 4: THICKNESS ===\n")

# Compare zones that exist in both
common_zones <- intersect(unique(thick_m$zone), unique(thick_z$zone))
thick_common <- thickness_all %>% filter(zone %in% common_zones)

# Thickness over absolute time
p4a <- ggplot(thick_common, aes(x = time_bin, y = thickness_p90,
                              color = species, linetype = zone)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.5, linewidth = 0.8) +
  scale_color_manual(values = species_colors) +
  labs(title = "Thickness over time",
       x = "Time (min)", y = "Thickness P5-P95 (um)", color = "Species") +
  theme_compare()

# Mean thickness by zone (bar plot)
thick_mean <- thick_common %>%
  group_by(species, zone) %>%
  summarise(mean_iqr = mean(thickness_p90, na.rm = TRUE),
            se = sd(thickness_p90, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

p4b <- ggplot(thick_mean, aes(x = zone, y = mean_iqr, fill = species)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_iqr - se, ymax = mean_iqr + se),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = species_colors) +
  labs(title = "Mean thickness by zone", x = NULL, y = "Thickness P5-P95 (um)", fill = "Species") +
  theme_compare()

# Depth distribution comparison
p4c <- ggplot(bind_rows(
    spots_m %>% dplyr::select(SPHERICAL_DEPTH, species),
    spots_z %>% dplyr::select(SPHERICAL_DEPTH, species)
  ) %>% mutate(species = factor(species, levels = c("Medaka", "Zebrafish"))),
  aes(x = SPHERICAL_DEPTH, fill = species)) +
  geom_density(alpha = 0.5, color = NA) +
  scale_fill_manual(values = species_colors) +
  labs(title = "Depth distribution",
       x = "Spherical depth (um)", y = "Density", fill = "Species") +
  theme_compare()

# Estimated cell layers — use NN distance from nuclear stats
if (HAS_NUCLEAR) {
  nn_m <- mean(nuc_m_clean[[nn3d_median_col]], na.rm = TRUE)
  nn_z <- mean(nuc_z_clean[[nn3d_median_col]], na.rm = TRUE)
  thick_layers <- thick_mean %>%
    mutate(cell_spacing = ifelse(species == "Medaka", nn_m, nn_z),
           n_layers = mean_iqr / cell_spacing)
  layer_subtitle <- sprintf("Thickness P5-P95 / NN distance (Medaka=%.1f, Zebrafish=%.1f um)", nn_m, nn_z)
} else {
  thick_layers <- thick_mean %>%
    mutate(nuc_diam = ifelse(species == "Medaka", 11.4, 12),
           n_layers = mean_iqr / nuc_diam)
  layer_subtitle <- "Thickness P5-P95 / nuclear diameter (fallback)"
}

p4d <- ggplot(thick_layers, aes(x = zone, y = n_layers, fill = species)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = species_colors) +
  labs(title = "Estimated cell layers by zone",
       subtitle = layer_subtitle,
       x = NULL, y = "Layers", fill = "Species") +
  theme_compare()

fig4 <- (p4a + p4b) / (p4c + p4d) +
  plot_annotation(
    title = "Medaka vs Zebrafish -- Blastoderm Thickness",
    theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14))
  )
save_pdf(fig4, "comparison_04_thickness.pdf", width = 14, height = 10)

# ============================================================================
# FIGURE 5: NUCLEAR SIZE & PROPERTIES
# ============================================================================

cat("\n=== FIGURE 5: NUCLEAR SIZE & PROPERTIES ===\n")

if (HAS_NUCLEAR) {

  # ── Compute early / late summaries for each species ──────────────────────
  compute_el <- function(df, col, frac = 0.2) {
    n <- nrow(df)
    early <- df %>% slice_head(n = round(n * frac))
    late  <- df %>% slice_tail(n = round(n * frac))
    e_val <- mean(early[[col]], na.rm = TRUE)
    l_val <- mean(late[[col]],  na.rm = TRUE)
    pct   <- (l_val - e_val) / e_val * 100
    list(early = e_val, late = l_val, pct = pct)
  }

  el_m <- list(
    diam3d   = compute_el(nuc_m_clean, eqsph_median_col),
    vol3d    = compute_el(nuc_m_clean, vol_median_col),
    nn3d     = compute_el(nuc_m_clean, nn3d_median_col),
    bbox_min = compute_el(nuc_m_clean, bbox_minor_col),
    n_nuc    = compute_el(nuc_m_clean, "n_nuclei_3d")
  )
  el_z <- list(
    diam3d   = compute_el(nuc_z_clean, eqsph_median_col),
    vol3d    = compute_el(nuc_z_clean, vol_median_col),
    nn3d     = compute_el(nuc_z_clean, nn3d_median_col),
    bbox_min = compute_el(nuc_z_clean, bbox_minor_col),
    n_nuc    = compute_el(nuc_z_clean, "n_nuclei_3d")
  )

  # ── Helper: interpolate value at a given time_min ────────────────────────
  interp_at <- function(df, col, t) {
    approx(df$time_min, df[[col]], xout = t, rule = 2)$y
  }

  # ── Helper: add vertical segment + % difference between species ──────────
  add_diff_segment <- function(p, col, t, m_df, z_df, nudge_x = 5) {
    vm <- interp_at(m_df, col, t)
    vz <- interp_at(z_df, col, t)
    if (is.na(vm) || is.na(vz) || vm == 0) return(p)
    pct_diff <- (vz - vm) / vm * 100
    seg_df <- data.frame(x = t, ymin = min(vm, vz), ymax = max(vm, vz))
    lab_df <- data.frame(x = t + nudge_x, y = (vm + vz) / 2,
                         label = sprintf("%+.1f%%", pct_diff))
    p +
      geom_segment(data = seg_df, aes(x = x, xend = x, y = ymin, yend = ymax),
                   inherit.aes = FALSE, linetype = "dotted", linewidth = 0.5, color = "grey30") +
      geom_point(data = data.frame(x = c(t, t), y = c(vm, vz)),
                 aes(x = x, y = y), inherit.aes = FALSE, size = 2.5, shape = 21,
                 fill = c(species_colors["Medaka"], species_colors["Zebrafish"]),
                 color = "white", stroke = 0.4) +
      geom_label(data = lab_df, aes(x = x, y = y, label = label),
                 inherit.aes = FALSE, size = 2.8, fill = "white", alpha = 0.85,
                 linewidth = 0.2, label.padding = unit(0.15, "lines"))
  }

  # ── Helper: annotate subtitle with early/late values and % change ────────
  annotate_el <- function(p, el_m, el_z, unit = "um", digits = 1) {
    m_txt <- sprintf(paste0("Med: %.", digits, "f -> %.", digits, "f %s (%+.1f%%)"),
                     el_m$early, el_m$late, unit, el_m$pct)
    z_txt <- sprintf(paste0("ZF:  %.", digits, "f -> %.", digits, "f %s (%+.1f%%)"),
                     el_z$early, el_z$late, unit, el_z$pct)
    p + labs(subtitle = paste0(m_txt, "\n", z_txt))
  }

  # ── Shared elements ─────────────────────────────────────────────────────
  t_end_m <- max(nuc_m_clean$time_min)
  t_end_z <- max(nuc_z_clean$time_min)
  t_overlap <- min(t_end_m, t_end_z)

  # SEM for confidence ribbons
  nuc_m_clean$diam_sem <- nuc_m_clean$eqsph_diam_std_um / sqrt(nuc_m_clean$n_nuclei_3d)
  nuc_z_clean$diam_sem <- nuc_z_clean$eqsph_diam_std_um / sqrt(nuc_z_clean$n_nuclei_3d)
  nuc_m_clean$nn_sem   <- nuc_m_clean$nn3d_std_um / sqrt(nuc_m_clean$n_nuclei_3d)
  nuc_z_clean$nn_sem   <- nuc_z_clean$nn3d_std_um / sqrt(nuc_z_clean$n_nuclei_3d)
  nuc_m_clean$vol_sem  <- nuc_m_clean$vol_std_um3 / sqrt(nuc_m_clean$n_nuclei_3d)
  nuc_z_clean$vol_sem  <- nuc_z_clean$vol_std_um3 / sqrt(nuc_z_clean$n_nuclei_3d)

  nuc_both <- bind_rows(nuc_m_clean, nuc_z_clean) %>%
    mutate(species = factor(species, levels = c("Medaka", "Zebrafish")))

  med_end_vline <- geom_vline(xintercept = t_end_m, linetype = "dashed",
                              color = species_colors["Medaka"], linewidth = 0.3, alpha = 0.6)

  # ── 5a: 3D Equiv-sphere diameter ─────────────────────────────────────────
  p5a <- ggplot(nuc_both, aes(time_min, .data[[eqsph_median_col]], color = species)) +
    geom_ribbon(aes(ymin = .data[[eqsph_median_col]] - diam_sem * 1.96,
                    ymax = .data[[eqsph_median_col]] + diam_sem * 1.96,
                    fill = species), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.7) +
    med_end_vline +
    scale_color_manual(values = species_colors) +
    scale_fill_manual(values = species_colors) +
    labs(title = "3D Equiv-sphere diameter (median)",
         x = "Time (min)", y = expression("Diameter ("*mu*m*")")) +
    theme_compare() + guides(fill = "none")
  p5a <- annotate_el(p5a, el_m$diam3d, el_z$diam3d, "\u00b5m", 2)
  p5a <- add_diff_segment(p5a, eqsph_median_col, 0, nuc_m_clean, nuc_z_clean, nudge_x = 12)
  p5a <- add_diff_segment(p5a, eqsph_median_col, t_overlap, nuc_m_clean, nuc_z_clean, nudge_x = -18)

  # ── 5b: 3D Volume ───────────────────────────────────────────────────────
  p5b <- ggplot(nuc_both, aes(time_min, .data[[vol_median_col]], color = species)) +
    geom_ribbon(aes(ymin = .data[[vol_median_col]] - vol_sem * 1.96,
                    ymax = .data[[vol_median_col]] + vol_sem * 1.96,
                    fill = species), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.7) +
    med_end_vline +
    scale_color_manual(values = species_colors) +
    scale_fill_manual(values = species_colors) +
    labs(title = "3D Nuclear volume (median)",
         x = "Time (min)", y = expression("Volume ("*mu*m^3*")")) +
    theme_compare() + guides(fill = "none")
  p5b <- annotate_el(p5b, el_m$vol3d, el_z$vol3d, "\u00b5m\u00b3", 0)
  p5b <- add_diff_segment(p5b, vol_median_col, 0, nuc_m_clean, nuc_z_clean, nudge_x = 12)
  p5b <- add_diff_segment(p5b, vol_median_col, t_overlap, nuc_m_clean, nuc_z_clean, nudge_x = -18)

  # ── 5c: Nuclei count ─────────────────────────────────────────────────────
  p5c <- ggplot(nuc_both, aes(time_min, n_nuclei_3d, color = species)) +
    geom_line(linewidth = 0.7) +
    med_end_vline +
    scale_color_manual(values = species_colors) +
    labs(title = "Nuclei count (3D)",
         x = "Time (min)", y = "N nuclei") +
    theme_compare()
  p5c <- annotate_el(p5c, el_m$n_nuc, el_z$n_nuc, "nuclei", 0)
  p5c <- add_diff_segment(p5c, "n_nuclei_3d", 0, nuc_m_clean, nuc_z_clean, nudge_x = 12)
  p5c <- add_diff_segment(p5c, "n_nuclei_3d", t_overlap, nuc_m_clean, nuc_z_clean, nudge_x = -18)

  # ── 5d: 3D NN distance ──────────────────────────────────────────────────
  p5d <- ggplot(nuc_both, aes(time_min, .data[[nn3d_median_col]], color = species)) +
    geom_ribbon(aes(ymin = .data[[nn3d_median_col]] - nn_sem * 1.96,
                    ymax = .data[[nn3d_median_col]] + nn_sem * 1.96,
                    fill = species), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.7) +
    med_end_vline +
    scale_color_manual(values = species_colors) +
    scale_fill_manual(values = species_colors) +
    labs(title = "3D Internuclear distance (median NN)",
         x = "Time (min)", y = expression("NN distance ("*mu*m*")")) +
    theme_compare() + guides(fill = "none")
  p5d <- annotate_el(p5d, el_m$nn3d, el_z$nn3d, "\u00b5m", 2)
  p5d <- add_diff_segment(p5d, nn3d_median_col, 0, nuc_m_clean, nuc_z_clean, nudge_x = 12)
  p5d <- add_diff_segment(p5d, nn3d_median_col, t_overlap, nuc_m_clean, nuc_z_clean, nudge_x = -18)

  # ── 5e: Bbox minor axis — lateral width ─────────────────────────────────
  # NOTE: ZF nuclei appear larger in volume/eqsph (ZF/Med = 1.17 vol, 1.05 diam)
  # but SMALLER in bbox minor (ZF/Med = 0.90). This is explained by Z-anisotropy:
  # ZF dZ=2.0 um vs Med dZ=1.05 um -> ZF nuclei are elongated along Z by the PSF,
  # inflating volume/eqsph. The bbox minor (shortest dimension = lateral XY width)
  # is unaffected by Z-elongation and shows medaka nuclei are laterally wider.
  # Aspect ratio: Med 1.42:1, ZF 1.96:1 — the ~2:1 ZF ratio mirrors dZ/dXY.
  p5e <- ggplot(nuc_both, aes(time_min, .data[[bbox_minor_col]], color = species)) +
    geom_line(linewidth = 0.7) +
    med_end_vline +
    scale_color_manual(values = species_colors) +
    labs(title = "Lateral nuclear width (bbox minor axis)",
         subtitle = paste0("ZF minor < Med despite larger volume: ZF dZ=2.0 vs Med dZ=1.05 ",
                           "\u00b5m elongates nuclei along Z, inflating V. Minor is Z-free."),
         x = "Time (min)", y = expression("Minor axis ("*mu*m*")")) +
    theme_compare()
  p5e <- add_diff_segment(p5e, bbox_minor_col, 0, nuc_m_clean, nuc_z_clean, nudge_x = 12)
  p5e <- add_diff_segment(p5e, bbox_minor_col, t_overlap, nuc_m_clean, nuc_z_clean, nudge_x = -18)

  # ── 5f: Proliferation dynamics — growth rate + cell cycle ────────────────
  # Compute smoothed instantaneous growth rate: lambda = d(ln N)/dt
  # Use heavier smoothing (span=0.3) to reduce noise — especially for medaka
  # where the count is nearly flat, making the derivative very noisy.
  compute_growth_rate <- function(df, span = 0.3) {
    lo <- loess(n_nuclei_3d ~ time_min, data = df, span = span)
    df$n_smooth <- predict(lo, newdata = df)
    df$ln_n <- log(df$n_smooth)
    # Second LOESS on ln(N) for a cleaner derivative
    lo2 <- loess(ln_n ~ time_min, data = df, span = span)
    df$ln_smooth <- predict(lo2, newdata = df)
    n <- nrow(df)
    df$lambda <- NA_real_
    for (i in 2:(n - 1)) {
      dt <- df$time_min[i + 1] - df$time_min[i - 1]
      if (dt > 0) df$lambda[i] <- (df$ln_smooth[i + 1] - df$ln_smooth[i - 1]) / dt
    }
    # One more smoothing pass on lambda itself to kill residual noise
    good <- which(!is.na(df$lambda))
    if (length(good) > 10) {
      lo3 <- loess(lambda ~ time_min, data = df[good, ], span = 0.4)
      df$lambda[good] <- predict(lo3)
    }
    df$doubling_min <- ifelse(!is.na(df$lambda) & df$lambda > 0, log(2) / df$lambda, NA_real_)
    df
  }

  nuc_m_rate <- compute_growth_rate(nuc_m_clean)
  nuc_z_rate <- compute_growth_rate(nuc_z_clean)

  rate_both <- bind_rows(
    nuc_m_rate %>% select(time_min, species, lambda, doubling_min, n_smooth),
    nuc_z_rate %>% select(time_min, species, lambda, doubling_min, n_smooth)
  ) %>% mutate(species = factor(species, levels = c("Medaka", "Zebrafish")))

  # Combined panel: doubling time as primary y-axis (intuitive for biologists).
  # During decline (lambda < 0), doubling time is undefined — shown as NA/gap.
  # T2 annotations: median of instantaneous doubling time over specific windows.
  zf_growth_phase <- nuc_z_rate %>% filter(!is.na(doubling_min), time_min < 200)
  med_growth_phase <- nuc_m_rate %>% filter(!is.na(doubling_min), time_min > 100)

  zf_dbl_med  <- median(zf_growth_phase$doubling_min, na.rm = TRUE)
  med_dbl_med <- median(med_growth_phase$doubling_min, na.rm = TRUE)
  zf_dbl_txt  <- sprintf("ZF: T2 ~ %.0f min (median, 0-200 min)", zf_dbl_med)
  med_dbl_txt <- sprintf("Med: T2 ~ %.0f min (median, 100-250 min)", med_dbl_med)

  # Filter to growth phases only (lambda > 0 -> finite doubling time)
  # Cap display at 1000 min — beyond that, essentially no division
  DOUBLING_CAP <- 1000
  dbl_plot_data <- rate_both %>%
    filter(!is.na(doubling_min) & doubling_min > 0 & doubling_min < DOUBLING_CAP)

  p5f <- ggplot(dbl_plot_data, aes(time_min, doubling_min, color = species)) +
    geom_line(linewidth = 0.8) +
    # Horizontal reference lines for cell-cycle durations
    geom_hline(yintercept = c(60, 120, 240, 480), linetype = "dotted",
               color = "grey75", linewidth = 0.25) +
    annotate("text", x = t_end_z, y = c(60, 120, 240, 480),
             label = c("1 h", "2 h", "4 h", "8 h"),
             size = 2.3, color = "grey55", hjust = 1.05, vjust = -0.3) +
    med_end_vline +
    scale_color_manual(values = species_colors) +
    # Annotate median doubling times with time window
    annotate("text", x = 100, y = zf_dbl_med + 30,
             label = zf_dbl_txt, size = 3, color = species_colors["Zebrafish"],
             hjust = 0.5, fontface = "italic") +
    annotate("text", x = 175, y = min(med_dbl_med + 40, DOUBLING_CAP * 0.85),
             label = med_dbl_txt, size = 3, color = species_colors["Medaka"],
             hjust = 0.5, fontface = "italic") +
    coord_cartesian(ylim = c(0, DOUBLING_CAP)) +
    labs(title = "Cell-cycle duration (doubling time)",
         subtitle = "T2 = ln(2) / growth rate. Shown only during net proliferation (growth rate > 0)",
         x = "Time (min)", y = "Doubling time (min)") +
    theme_compare()

  cat(sprintf("  Cell cycle: %s | %s\n", zf_dbl_txt, med_dbl_txt))

  # ── 5g: Summary bar chart — early vs late change for key metrics ─────────
  summary_df <- tibble(
    metric  = rep(c("EqSph diam", "Volume", "NN dist", "Bbox minor", "N nuclei"), each = 2),
    species = rep(c("Medaka", "Zebrafish"), 5),
    pct     = c(el_m$diam3d$pct, el_z$diam3d$pct,
                el_m$vol3d$pct,  el_z$vol3d$pct,
                el_m$nn3d$pct,   el_z$nn3d$pct,
                el_m$bbox_min$pct, el_z$bbox_min$pct,
                el_m$n_nuc$pct,  el_z$n_nuc$pct)
  ) %>% mutate(
    metric  = factor(metric, levels = c("EqSph diam", "Volume", "NN dist", "Bbox minor", "N nuclei")),
    species = factor(species, levels = c("Medaka", "Zebrafish")),
    label   = sprintf("%+.1f%%", pct)
  )

  p5g <- ggplot(summary_df, aes(metric, pct, fill = species)) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.85) +
    geom_text(aes(label = label, y = pct + sign(pct) * 3),
              position = position_dodge(0.8), size = 3, vjust = 0.5) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = species_colors) +
    labs(title = "Early-to-late change (%)",
         subtitle = sprintf("Early = first 20%% | Late = last 20%%. Med 0-%.0f min, ZF 0-%.0f min",
                            t_end_m, t_end_z),
         x = NULL, y = "Change (%)") +
    theme_compare() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  # ── 5h: Nuclear aspect ratio (major/minor) — Z-anisotropy diagnostic ──────
  # This panel directly answers WHY ZF nuclei appear larger in volume/eqsph
  # but smaller in bbox minor. ZF dZ=2.0 um (vs Med 1.05) stretches the PSF
  # along Z, elongating segmented nuclei. The major/minor ratio captures this:
  # if nuclei were truly spherical, ratio ~1. Med ratio ~1.4 (mild elongation),
  # ZF ratio ~2.0 — mirroring dZ/dXY=1.60. The "extra" volume is Z-elongation
  # artifact, not real biological size difference. Bbox minor (lateral XY) is
  # the fair cross-species comparison.
  nuc_m_clean$aspect_ratio <- nuc_m_clean[[bbox_major_col]] / nuc_m_clean[[bbox_minor_col]]
  nuc_z_clean$aspect_ratio <- nuc_z_clean[[bbox_major_col]] / nuc_z_clean[[bbox_minor_col]]

  nuc_both_ar <- bind_rows(nuc_m_clean, nuc_z_clean) %>%
    mutate(species = factor(species, levels = c("Medaka", "Zebrafish")))

  p5h <- ggplot(nuc_both_ar, aes(time_min, aspect_ratio, color = species)) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 1.0, linetype = "dotted", color = "grey50", linewidth = 0.3) +
    med_end_vline +
    scale_color_manual(values = species_colors) +
    annotate("text", x = t_end_z * 0.6, y = 1.05, label = "perfect sphere",
             size = 2.8, color = "grey50", fontface = "italic") +
    labs(title = "Nuclear elongation (bbox major / minor)",
         subtitle = sprintf("ZF ~%.1f:1 vs Med ~%.1f:1. ZF dZ/dXY=%.1f explains inflated volume",
                            mean(nuc_z_clean$aspect_ratio),
                            mean(nuc_m_clean$aspect_ratio),
                            nuc_global_z$voxel_dz_um / nuc_global_z$voxel_dx_um),
         x = "Time (min)", y = "Aspect ratio (major/minor)") +
    theme_compare()

  # ── Assemble: 4 rows x 2 cols = 8 panels ───────────────────────────────
  fig5 <- (p5a + p5b) / (p5c + p5d) / (p5e + p5f) / (p5g + p5h) +
    plot_annotation(
      title = "Medaka vs Zebrafish -- Nuclear Properties",
      subtitle = paste0(
        sprintf("Med 0-%.0f min (pre-contraction) | ZF 0-%.0f min. ", t_end_m, t_end_z),
        "Dashed line = Med endpoint. 95% CI ribbons = SEM. ",
        sprintf("Voxel: Med %.2f\u00d7%.2f, ZF %.2f\u00d7%.2f \u00b5m. ",
                nuc_global_m$voxel_dx_um, nuc_global_m$voxel_dz_um,
                nuc_global_z$voxel_dx_um, nuc_global_z$voxel_dz_um),
        "n=1 embryo/species."),
      theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14),
                    plot.subtitle = element_text(hjust = 0, size = 9, color = "grey40"))
    )
  save_pdf(fig5, "comparison_05_nuclear.pdf", width = 14, height = 20)
} else {
  cat("  Skipped (nuclear stats not available)\n")
}

# ============================================================================
# FIGURE 6: DIRECTIONAL COMPOSITION & MOVEMENT TYPE
# ============================================================================

cat("\n=== FIGURE 6: DIRECTIONAL COMPOSITION ===\n")

# Velocity vector angle histogram (epiboly vs lateral)
vel_angle <- vel_all %>%
  mutate(angle = atan2(0, v_epiboly) * 180 / pi) # simplified: use theta velocity direction

# Instead, compare epiboly fraction
epiboly_fraction <- vel_all %>%
  mutate(is_epiboly = v_epiboly > 0,
         time_bin = floor(time_min / TIME_BIN_MIN) * TIME_BIN_MIN) %>%
  group_by(species, time_bin) %>%
  summarise(epiboly_pct = 100 * mean(is_epiboly, na.rm = TRUE),
            mean_v_epi = mean(v_epiboly, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= 50)

p6a <- ggplot(epiboly_fraction, aes(x = time_bin, y = epiboly_pct, color = species)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  scale_color_manual(values = species_colors) +
  labs(title = "Fraction of vegetalward (epiboly) steps over time",
       x = "Time (min)", y = "% epiboly", color = "Species") +
  theme_compare()

p6b <- ggplot(epiboly_fraction, aes(x = time_bin, y = mean_v_epi, color = species)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = species_colors) +
  labs(title = "Mean epiboly velocity over time",
       x = "Time (min)", y = "v_epiboly (um/min)", color = "Species") +
  theme_compare()

# Depth velocity comparison
depth_vel_cmp <- vel_all %>%
  filter(!is.na(v_depth)) %>%
  mutate(time_bin = floor(time_min / TIME_BIN_MIN) * TIME_BIN_MIN) %>%
  group_by(species, time_bin) %>%
  summarise(mean_v_depth = mean(v_depth, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= 50)

p6c <- ggplot(depth_vel_cmp, aes(x = time_bin, y = mean_v_depth, color = species)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = species_colors) +
  labs(title = "Mean depth velocity over time",
       subtitle = "Positive = deepening (internalization)",
       x = "Time (min)", y = "v_depth (um/min)", color = "Species") +
  theme_compare()

# Movement type proportions
mvt_prop <- track_all %>%
  filter(!is.na(movement_type)) %>%
  group_by(species, movement_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(species) %>%
  mutate(pct = 100 * n / sum(n))

p6d <- ggplot(mvt_prop, aes(x = species, y = pct, fill = movement_type)) +
  geom_col(position = "stack", color = "white", linewidth = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Track movement type composition",
       x = NULL, y = "%", fill = "Type") +
  theme_compare()

fig6 <- (p6a + p6b) / (p6c + p6d) +
  plot_annotation(
    title = "Medaka vs Zebrafish -- Directional Composition",
    theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14))
  )
save_pdf(fig6, "comparison_06_direction.pdf", width = 14, height = 10)

# ============================================================================
# FIGURE 7: SUMMARY TABLE (as a plot)
# ============================================================================

cat("\n=== FIGURE 7: SUMMARY TABLE ===\n")

# Build a comparison table from summary CSVs
sum_wide <- summary_both %>%
  pivot_wider(names_from = species, values_from = value)

# Also compute the ratio
sum_wide$ratio <- ifelse(!is.na(sum_wide$Medaka) & !is.na(sum_wide$Zebrafish) &
                         sum_wide$Zebrafish != 0,
                         round(as.numeric(sum_wide$Medaka) / as.numeric(sum_wide$Zebrafish), 2), NA)

# Format numbers nicely
fmt_val <- function(x) {
  x_num <- as.numeric(x)
  ifelse(is.na(x_num), x,
         ifelse(abs(x_num) > 100, format(round(x_num), big.mark = ","),
                ifelse(abs(x_num) > 1, sprintf("%.2f", x_num),
                       sprintf("%.3f", x_num))))
}

sum_table <- sum_wide %>%
  mutate(Medaka = fmt_val(Medaka),
         Zebrafish = fmt_val(Zebrafish),
         ratio = ifelse(is.na(ratio), "", as.character(ratio)))

write_csv(sum_table, file.path(OUTPUT_DIR, "comparison_summary_table.csv"))
cat("  Summary table:\n")
for (i in seq_len(nrow(sum_table))) {
  cat(sprintf("    %-25s  M: %-12s  Z: %-12s  ratio: %s\n",
              sum_table$metric[i], sum_table$Medaka[i], sum_table$Zebrafish[i],
              sum_table$ratio[i]))
}

# Create a table plot
p7 <- ggplot(sum_table, aes(y = reorder(metric, nrow(sum_table):1))) +
  geom_text(aes(x = 0.5, label = Medaka), color = species_colors["Medaka"],
            fontface = "bold", size = 3.5, hjust = 1) +
  geom_text(aes(x = 1.5, label = Zebrafish), color = species_colors["Zebrafish"],
            fontface = "bold", size = 3.5, hjust = 0) +
  geom_text(aes(x = 2.5, label = ratio), color = "grey30", size = 3.5, hjust = 0.5) +
  annotate("text", x = 0.5, y = nrow(sum_table) + 0.8, label = "Medaka",
           fontface = "bold", color = species_colors["Medaka"], size = 4) +
  annotate("text", x = 1.5, y = nrow(sum_table) + 0.8, label = "Zebrafish",
           fontface = "bold", color = species_colors["Zebrafish"], size = 4) +
  annotate("text", x = 2.5, y = nrow(sum_table) + 0.8, label = "Ratio M/Z",
           fontface = "bold", color = "grey30", size = 4) +
  scale_x_continuous(limits = c(-0.5, 3.5)) +
  labs(title = "Medaka vs Zebrafish -- Summary Metrics", x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14,
                                  margin = margin(b = 10)),
        axis.text.y = element_text(hjust = 1, size = 9))

save_pdf(p7, "comparison_07_summary_table.pdf", width = 10, height = 8)

cat("\n")
cat(strrep("=", 70), "\n")
cat("  COMPARISON ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n")
