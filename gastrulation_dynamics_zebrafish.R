# =============================================================================
# Gastrulation Dynamics — ZEBRAFISH Focused Analysis
# =============================================================================
#
# Same classification algorithm as medaka. Key biological differences:
#   - Epiboly-dominated movement (convergence naturally minimal/absent)
#   - Frame interval = 120 sec (2 min)
#
# Core analyses: Velocity/Acceleration, Speed, Straightness, Persistence,
# Neighbor coordination, Directional decomposition, MSD, Thickness
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

SPECIES <- "Zebrafish"
FRAME_INTERVAL_SEC <- 120
FRAME_INTERVAL_MIN <- FRAME_INTERVAL_SEC / 60

INPUT_DIR  <- "oriented_zebrafish_ultrack"
OUTPUT_DIR <- "analysis_output_zebrafish"

MSD_MAX_LAG    <- 30     # fewer lags — frames are 2 min apart
MSD_N_TRACKS   <- 10000
NEIGHBOR_MAX_DIST <- 100
NEIGHBOR_DIST_BIN <- 10
THETA_BIN      <- 3
PHI_BIN        <- 5
FLOW_MIN_N     <- 10
TILE_DTHETA    <- 2.0
TILE_DPHI      <- 3.0
MIN_TRACK_LENGTH <- 5
DEPTH_EPOCH_SIZE <- 10   # frames (~20 min)
INGRESSION_DEPTH_PCTILE <- 0.95
INGRESSION_ONSET_GAP_WINDOW <- 5
ZONE_MARGIN_WIDTH  <- 8  # wider zones for zebrafish (larger theta range)
ZONE_PREMARG_WIDTH <- 15
SMOOTH_K       <- 3    # position smoothing window (frames); 3×2min = 6 min (comparable to medaka 5×0.5min)
VELOCITY_LAG   <- 1    # frames over which to compute velocity (1×120s = 2 min)
NUCLEAR_DIR    <- "nuclei_stats_zebrafish"

# Voxel calibration — ultrack outputs positions in pixels (voxel indices).
# Voxel sizes from TIFF metadata (nuclear_stats_global.tsv): dXY=1.24785 um, dZ=1.25 um.
# Near-isotropic (<0.2% anisotropy), so a single scale factor is applied after rotation.
VOXEL_SIZE_UM  <- 1.24785  # um per pixel (XY)

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# =============================================================================
# THEME & HELPERS
# =============================================================================

theme_dynamics <- function(base_size = 10) {
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

theme_spatial <- function(base_size = 10) {
  theme_dynamics(base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "grey15", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
}

save_pdf <- function(plot, filename, width = 12, height = 8) {
  path <- file.path(OUTPUT_DIR, filename)
  ggsave(path, plot, width = width, height = height, device = "pdf")
  cat(sprintf("  -> Saved %s\n", path))
}

frame_to_min <- function(f) f * FRAME_INTERVAL_MIN

zone_colors <- c(
  "Animal cap"       = "#2166AC",
  "Pre-marginal"     = "#67A9CF",
  "Margin layer"     = "#EF8A62",
  "Sub-marginal"     = "#B2182B"
)

ZONE_LEVELS <- c("Animal cap", "Pre-marginal", "Margin layer", "Sub-marginal")

# Zebrafish step types: same as movement types (4 categories, matching medaka)
step_colors <- c(
  "Epiboly"      = "#EF8A62",
  "Animalward"   = "#FDDBC7",
  "Convergence"  = "#1B7837"
)

mvt_colors <- c(
  "Epiboly"      = "#EF8A62",
  "Animalward"   = "#FDDBC7",
  "Convergence"  = "#1B7837",
  "Ingression"   = "#2166AC"
)

THETA_LABEL <- "Theta (deg): 0=animal pole, 90=equator"
PHI_LABEL   <- "Phi (deg)"

# ###########################################################################
# PHASE 1: DATA LOADING & COMPUTATIONS
# ###########################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("  GASTRULATION DYNAMICS — ZEBRAFISH\n")
cat(strrep("=", 70), "\n")

# =============================================================================
# STEP 1: LOAD DATA
# =============================================================================

cat("\n=== STEP 1: LOADING DATA ===\n")

spots <- fread(file.path(INPUT_DIR, "oriented_tracks_zebrafish.csv"), showProgress = FALSE)
setDF(spots)

sphere_raw <- read_csv(file.path(INPUT_DIR, "sphere_params.csv"), show_col_types = FALSE)
sphere <- setNames(sphere_raw$value, sphere_raw$parameter)
SPHERE_RADIUS <- as.numeric(sphere["radius"])

# --- Voxel calibration: convert pixel positions to microns ---
cat(sprintf("  Applying voxel calibration: %.5f um/px\n", VOXEL_SIZE_UM))
for (col in c("POSITION_X", "POSITION_Y", "POSITION_Z", "RADIAL_DIST", "SPHERICAL_DEPTH")) {
  spots[[col]] <- spots[[col]] * VOXEL_SIZE_UM
}
SPHERE_RADIUS <- SPHERE_RADIUS * VOXEL_SIZE_UM

cat(sprintf("  Spots: %s rows | Tracks: %s | Frames: %d-%d\n",
            format(nrow(spots), big.mark = ","),
            format(length(unique(spots$TRACK_ID)), big.mark = ","),
            min(spots$FRAME), max(spots$FRAME)))
cat(sprintf("  Sphere radius: %.1f um\n", SPHERE_RADIUS))
cat(sprintf("  Theta range: %.1f - %.1f | Phi range: %.1f - %.1f\n",
            min(spots$THETA_DEG), max(spots$THETA_DEG),
            min(spots$PHI_DEG), max(spots$PHI_DEG)))

# Load landmarks
LANDMARK_FILE <- file.path(INPUT_DIR, "gastrulation_landmarks.csv")
MARGIN_THETA_LANDMARK <- NA_real_
INGRESSION_THETA_LM   <- NA_real_
INGRESSION_PHI_LM     <- NA_real_

if (file.exists(LANDMARK_FILE)) {
  lm_df <- read_csv(LANDMARK_FILE, show_col_types = FALSE)
  margin_row <- lm_df %>% filter(landmark == "margin")
  if (nrow(margin_row) == 1 && "theta" %in% names(margin_row)) {
    MARGIN_THETA_LANDMARK <- margin_row$theta[1]
    cat(sprintf("  Landmark margin theta: %.1f\n", MARGIN_THETA_LANDMARK))
  }
  ingr_row <- lm_df %>% filter(landmark == "ingression_center")
  if (nrow(ingr_row) == 1 && "theta" %in% names(ingr_row)) {
    INGRESSION_THETA_LM <- ingr_row$theta[1]
    INGRESSION_PHI_LM   <- ingr_row$phi[1]
    cat(sprintf("  Landmark ingression: theta=%.1f, phi=%.1f\n",
                INGRESSION_THETA_LM, INGRESSION_PHI_LM))
  }
}

spots <- spots %>% mutate(time_min = frame_to_min(FRAME))

# =============================================================================
# STEP 2: MARGIN DETECTION
# =============================================================================

cat("\n=== STEP 2: MARGIN DETECTION ===\n")

# Margin: use landmark if available, else auto-detect
if (!is.na(MARGIN_THETA_LANDMARK)) {
  MARGIN_THETA <- MARGIN_THETA_LANDMARK
  cat(sprintf("  Using landmark margin: %.1f deg\n", MARGIN_THETA))
} else {
  MARGIN_THETA <- quantile(spots$THETA_DEG, 0.90)
  cat(sprintf("  Auto-detected margin (90th pctile): %.1f deg\n", MARGIN_THETA))
}

# Dynamic margin tracking
WIN <- 10
margin_raw <- lapply(seq(0, max(spots$FRAME), WIN), function(t_start) {
  t_end <- min(t_start + WIN - 1, max(spots$FRAME))
  sub <- spots %>% filter(FRAME >= t_start & FRAME <= t_end)
  if (nrow(sub) < 100) return(NULL)
  dens <- sub %>%
    mutate(theta_bin = floor(THETA_DEG) + 0.5) %>%
    group_by(theta_bin) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(theta_bin) %>%
    mutate(n_smooth = zoo::rollmean(n, k = 5, fill = NA, align = "center"))
  dens$n_smooth[is.na(dens$n_smooth)] <- dens$n[is.na(dens$n_smooth)]
  peak_n <- max(dens$n_smooth, na.rm = TRUE)
  margin_dens <- max(dens$theta_bin[dens$n_smooth >= peak_n * 0.20])
  data.frame(FRAME_mid = (t_start + t_end) / 2,
             time_min  = frame_to_min((t_start + t_end) / 2),
             margin_raw = margin_dens)
}) %>% bind_rows()

if (nrow(margin_raw) >= 5) {
  lo <- loess(margin_raw ~ FRAME_mid, data = margin_raw, span = 0.75)
  margin_raw$margin_loess <- predict(lo)
} else {
  margin_raw$margin_loess <- margin_raw$margin_raw
}
margin_raw$margin_mono <- isoreg(margin_raw$margin_loess)$yf

margin_per_frame <- data.frame(FRAME = 0:max(spots$FRAME)) %>%
  mutate(time_min = frame_to_min(FRAME))
margin_interp <- approx(x = margin_raw$FRAME_mid, y = margin_raw$margin_mono,
                        xout = margin_per_frame$FRAME, rule = 2)
margin_per_frame$margin_smooth <- margin_interp$y
margin_per_frame$dmargin_dt <- c(NA, diff(margin_per_frame$margin_smooth)) / FRAME_INTERVAL_MIN

cat(sprintf("  Dynamic margin: %.1f -> %.1f deg\n",
            margin_per_frame$margin_smooth[1], tail(margin_per_frame$margin_smooth, 1)))

# =============================================================================
# STEP 2b: INGRESSION CLUSTER DETECTION (same algorithm as medaka)
# =============================================================================

cat("\n=== STEP 2b: INGRESSION CLUSTER DETECTION ===\n")

track_depth_all <- spots %>%
  group_by(TRACK_ID) %>%
  summarise(max_depth = max(SPHERICAL_DEPTH, na.rm = TRUE),
            mean_theta = mean(THETA_DEG, na.rm = TRUE),
            mean_phi   = mean(PHI_DEG, na.rm = TRUE),
            n_frames = n(), .groups = "drop") %>%
  filter(n_frames >= 3)

depth_threshold <- quantile(track_depth_all$max_depth, INGRESSION_DEPTH_PCTILE)
ingr_candidates <- track_depth_all %>% filter(max_depth >= depth_threshold)

if (nrow(ingr_candidates) >= 10) {
  INGRESSION_CENTER_THETA <- mean(ingr_candidates$mean_theta)
  INGRESSION_CENTER_PHI   <- mean(ingr_candidates$mean_phi)
  ingr_candidates <- ingr_candidates %>%
    mutate(ang_dist = sqrt((mean_theta - INGRESSION_CENTER_THETA)^2 +
                           (mean_phi   - INGRESSION_CENTER_PHI)^2))
  INGRESSION_RADIUS <- min(2 * sd(ingr_candidates$ang_dist) + 2, 15)
  INGRESSION_DEPTH_THRESH <- depth_threshold
} else {
  INGRESSION_CENTER_THETA <- ifelse(!is.na(INGRESSION_THETA_LM), INGRESSION_THETA_LM, 111)
  INGRESSION_CENTER_PHI   <- ifelse(!is.na(INGRESSION_PHI_LM), INGRESSION_PHI_LM, 85)
  INGRESSION_RADIUS <- 10
  INGRESSION_DEPTH_THRESH <- 60
}

cat(sprintf("  Cluster: theta=%.1f, phi=%.1f, radius=%.1f deg\n",
            INGRESSION_CENTER_THETA, INGRESSION_CENTER_PHI, INGRESSION_RADIUS))

# Onset detection from depth gap
depth_gap_series <- spots %>%
  mutate(dist_to_ingr_tmp = sqrt((THETA_DEG - INGRESSION_CENTER_THETA)^2 +
                                  (PHI_DEG - INGRESSION_CENTER_PHI)^2),
         in_zone_tmp = dist_to_ingr_tmp < INGRESSION_RADIUS,
         epoch = floor(FRAME / 10) * 10) %>%
  group_by(epoch, in_zone_tmp) %>%
  summarise(mean_depth = mean(SPHERICAL_DEPTH, na.rm = TRUE), n = n(), .groups = "drop") %>%
  pivot_wider(names_from = in_zone_tmp, values_from = c(mean_depth, n)) %>%
  mutate(gap = mean_depth_TRUE - mean_depth_FALSE) %>%
  filter(!is.na(gap)) %>% arrange(epoch)

gap_vec <- depth_gap_series$gap
n_pts <- length(gap_vec)
w <- INGRESSION_ONSET_GAP_WINDOW
slopes <- rep(NA_real_, n_pts)
for (i in (w + 1):(n_pts - w)) {
  idx <- (i - w):(i + w)
  slopes[i] <- coef(lm(gap_vec[idx] ~ seq_along(idx)))[2]
}

positive_slopes <- which(slopes > 0 & !is.na(slopes))
onset_idx <- NA_integer_
if (length(positive_slopes) >= 3) {
  for (j in seq_len(length(positive_slopes) - 2)) {
    if (positive_slopes[j + 1] == positive_slopes[j] + 1 &&
        positive_slopes[j + 2] == positive_slopes[j] + 2) {
      onset_idx <- positive_slopes[j]; break
    }
  }
}

if (!is.na(onset_idx)) {
  INGRESSION_ONSET_FRAME <- depth_gap_series$epoch[onset_idx]
} else {
  baseline_gap <- mean(gap_vec[1:min(5, n_pts)])
  baseline_sd  <- sd(gap_vec[1:min(10, n_pts)])
  sig <- which(gap_vec > baseline_gap + 2 * baseline_sd)
  INGRESSION_ONSET_FRAME <- if (length(sig) > 0) depth_gap_series$epoch[sig[1]] else round(max(spots$FRAME) * 0.3)
}

cat(sprintf("  Onset frame: %d (%.1f min)\n", INGRESSION_ONSET_FRAME,
            INGRESSION_ONSET_FRAME * FRAME_INTERVAL_MIN))

track_depth_post <- spots %>%
  filter(FRAME >= INGRESSION_ONSET_FRAME) %>%
  group_by(TRACK_ID) %>%
  summarise(max_depth_post = max(SPHERICAL_DEPTH, na.rm = TRUE),
            n_frames_post = n(),
            depth_change_post = last(SPHERICAL_DEPTH) - first(SPHERICAL_DEPTH),
            .groups = "drop") %>%
  filter(n_frames_post >= 3)

# =============================================================================
# STEP 3: VELOCITY & ACCELERATION
# =============================================================================

cat("\n=== STEP 3: COMPUTING VELOCITY & ACCELERATION ===\n")
cat(sprintf("  Smoothing: k=%d | Velocity lag: %d frames (%.0f s = %.1f min)\n",
            SMOOTH_K, VELOCITY_LAG, VELOCITY_LAG * FRAME_INTERVAL_SEC,
            VELOCITY_LAG * FRAME_INTERVAL_MIN))

# Use data.table for speed
setDT(spots)
setkey(spots, TRACK_ID, FRAME)

# Rolling mean smoothing using data.table's frollmean (C-level, fast)
for (col in c("POSITION_X", "POSITION_Y", "POSITION_Z", "RADIAL_DIST", "THETA_DEG", "PHI_DEG")) {
  sm_col <- paste0(col, "_SM")
  spots[, (sm_col) := frollmean(get(col), n = SMOOTH_K, align = "center"), by = TRACK_ID]
  spots[is.na(get(sm_col)), (sm_col) := get(col)]
}

# Velocity computation — shift by VELOCITY_LAG frames
VL <- as.integer(VELOCITY_LAG)
spots[, `:=`(
  dx        = POSITION_X_SM - shift(POSITION_X_SM, VL),
  dy        = POSITION_Y_SM - shift(POSITION_Y_SM, VL),
  dz        = POSITION_Z_SM - shift(POSITION_Z_SM, VL),
  dt_frames = FRAME - shift(FRAME, VL),
  dR        = RADIAL_DIST_SM - shift(RADIAL_DIST_SM, VL),
  dDepth    = SPHERICAL_DEPTH - shift(SPHERICAL_DEPTH, VL),
  dTheta    = THETA_DEG_SM - shift(THETA_DEG_SM, VL),
  dPhi      = PHI_DEG_SM - shift(PHI_DEG_SM, VL)
), by = TRACK_ID]

spots[, dPhi := fifelse(abs(dPhi) > 180, dPhi - sign(dPhi) * 360, dPhi)]
spots[, `:=`(
  disp_3d    = sqrt(dx^2 + dy^2 + dz^2),
  theta_rad  = THETA_DEG * pi / 180,
  dTheta_rad = dTheta * pi / 180,
  dPhi_rad   = dPhi * pi / 180
)]
spots[, `:=`(
  inst_speed     = disp_3d / (dt_frames * FRAME_INTERVAL_MIN),
  v_radial       = dR / (dt_frames * FRAME_INTERVAL_MIN),
  v_depth        = dDepth / (dt_frames * FRAME_INTERVAL_MIN),
  v_theta        = dTheta / (dt_frames * FRAME_INTERVAL_MIN),
  v_phi          = dPhi / (dt_frames * FRAME_INTERVAL_MIN),
  v_tangential   = RADIAL_DIST * sqrt(dTheta_rad^2 + (sin(theta_rad) * dPhi_rad)^2) /
                   (dt_frames * FRAME_INTERVAL_MIN),
  # Epiboly = vegetalward along meridian
  v_epiboly      = RADIAL_DIST * dTheta_rad / (dt_frames * FRAME_INTERVAL_MIN),
  # Convergence = azimuthal movement (same decomposition as medaka)
  v_convergence  = RADIAL_DIST * sin(theta_rad) * dPhi_rad /
                   (dt_frames * FRAME_INTERVAL_MIN)
)]

# Acceleration & turning angle (need lag of computed columns — use VL step spacing)
spots[, `:=`(
  d_speed  = inst_speed - shift(inst_speed, VL),
  dt2      = dt_frames + shift(dt_frames, VL),
  lag_dx   = shift(dx, VL),
  lag_dy   = shift(dy, VL),
  lag_dz   = shift(dz, VL),
  mag_prev = shift(disp_3d, VL)
), by = TRACK_ID]

spots[, `:=`(
  acceleration  = 2 * d_speed / (dt2 * FRAME_INTERVAL_MIN),
  dot_prev      = dx * lag_dx + dy * lag_dy + dz * lag_dz,
  mag_curr      = disp_3d
)]
spots[, cos_angle := dot_prev / (mag_curr * mag_prev)]
spots[, turning_angle := acos(pmin(pmax(cos_angle, -1), 1)) * 180 / pi]

# Clean up temp columns and convert back to data.frame
spots[, c("lag_dx", "lag_dy", "lag_dz", "d_speed", "dt2", "dot_prev",
          "mag_curr", "mag_prev", "cos_angle") := NULL]
setDF(spots)

spots_vel <- spots %>% filter(!is.na(dt_frames) & dt_frames == VELOCITY_LAG)

cat(sprintf("  Velocity steps: %s\n", format(nrow(spots_vel), big.mark = ",")))

# =============================================================================
# STEP 3b: SMOOTHING SENSITIVITY — TURNING ANGLE ARTIFACT TEST
# =============================================================================
# Test whether turning angles are artifacts of segmentation noise ("wiggling").
# If aggressive trajectory smoothing substantially reduces turning angles, the
# original values reflect noise.  If they persist, they reflect biological
# behaviour.  Inspired by wavelet-denoising approaches for trajectory
# smoothing (cf. DWT sym4 single-level approximation reconstruction).
#
# Method: re-compute turning angles from positions smoothed at four levels
# (raw → current → moderate → aggressive rolling mean), then compare the
# distributions.  Speed is monitored as a control and should remain stable.
# =============================================================================

cat("\n=== STEP 3b: SMOOTHING SENSITIVITY — TURNING ANGLE ARTIFACT TEST ===\n")

SMOOTH_LEVELS <- c(1L, as.integer(SMOOTH_K),
                    as.integer(2 * SMOOTH_K + 1),
                    as.integer(4 * SMOOTH_K + 1))
SMOOTH_LABELS <- c("Raw (k=1)",
                    sprintf("Current (k=%d)", SMOOTH_K),
                    sprintf("Moderate (k=%d)", 2 * SMOOTH_K + 1),
                    sprintf("Aggressive (k=%d)", 4 * SMOOTH_K + 1))

setDT(spots)
setkey(spots, TRACK_ID, FRAME)

smooth_val_results <- lapply(seq_along(SMOOTH_LEVELS), function(lvl) {
  k <- SMOOTH_LEVELS[lvl]
  label <- SMOOTH_LABELS[lvl]

  tmp <- spots[, .(TRACK_ID, FRAME, POSITION_X, POSITION_Y, POSITION_Z)]
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

  VL_int <- as.integer(VELOCITY_LAG)
  tmp[, `:=`(dx = sx - shift(sx, VL_int),
             dy = sy - shift(sy, VL_int),
             dz = sz - shift(sz, VL_int),
             dt_f = FRAME - shift(FRAME, VL_int)), by = TRACK_ID]
  tmp[, disp_3d := sqrt(dx^2 + dy^2 + dz^2)]
  tmp[, inst_speed := disp_3d / (dt_f * FRAME_INTERVAL_MIN)]

  tmp[, `:=`(lag_dx = shift(dx, VL_int),
             lag_dy = shift(dy, VL_int),
             lag_dz = shift(dz, VL_int),
             mag_prev = shift(disp_3d, VL_int)), by = TRACK_ID]
  tmp[, dot_prev := dx * lag_dx + dy * lag_dy + dz * lag_dz]
  tmp[, cos_ang := dot_prev / (disp_3d * mag_prev)]
  tmp[, turn_ang := acos(pmin(pmax(cos_ang, -1), 1)) * 180 / pi]

  valid <- tmp[!is.na(dt_f) & dt_f == VL_int & !is.na(turn_ang) & !is.na(inst_speed)]
  data.table(smooth_level = label, smooth_k = k,
             turning_angle = valid$turn_ang, inst_speed = valid$inst_speed)
})
smooth_val_results <- rbindlist(smooth_val_results)
smooth_val_results[, smooth_level := factor(smooth_level, levels = SMOOTH_LABELS)]

setDF(spots)

smooth_val_summary <- smooth_val_results[, .(
  median_turn  = median(turning_angle, na.rm = TRUE),
  mean_turn    = mean(turning_angle, na.rm = TRUE),
  sd_turn      = sd(turning_angle, na.rm = TRUE),
  median_speed = median(inst_speed, na.rm = TRUE),
  mean_speed   = mean(inst_speed, na.rm = TRUE),
  n = .N
), by = .(smooth_level, smooth_k)]

raw_median_turn  <- smooth_val_summary[smooth_k == 1, median_turn]
raw_median_speed <- smooth_val_summary[smooth_k == 1, median_speed]
smooth_val_summary[, pct_change_turn  := (median_turn  - raw_median_turn)  / raw_median_turn  * 100]
smooth_val_summary[, pct_change_speed := (median_speed - raw_median_speed) / raw_median_speed * 100]

cat("  Smoothing level sensitivity:\n")
for (i in seq_len(nrow(smooth_val_summary))) {
  s <- smooth_val_summary[i]
  cat(sprintf("    %s: median turn = %.1f deg (%+.1f%%), median speed = %.3f um/min (%+.1f%%) [n=%s]\n",
              s$smooth_level, s$median_turn, s$pct_change_turn,
              s$median_speed, s$pct_change_speed,
              format(s$n, big.mark = ",")))
}

aggr_turn_pct  <- smooth_val_summary[smooth_k == max(SMOOTH_LEVELS), pct_change_turn]
aggr_speed_pct <- smooth_val_summary[smooth_k == max(SMOOTH_LEVELS), pct_change_speed]
cat(sprintf("  CONCLUSION: Turning angle change (raw → aggressive): %+.1f%%\n", aggr_turn_pct))
cat(sprintf("              Speed change (raw → aggressive):         %+.1f%%\n", aggr_speed_pct))
if (abs(aggr_turn_pct) < 20) {
  cat("  → Turning angles are ROBUST to smoothing — unlikely to be noise artifacts.\n")
} else {
  cat("  → Turning angles are SENSITIVE to smoothing — may partially reflect segmentation noise.\n")
}

# =============================================================================
# STEP 4: ZONE ASSIGNMENT
# =============================================================================

cat("\n=== STEP 4: ZONE ASSIGNMENT ===\n")

ZONE_MARGIN_TOP  <- MARGIN_THETA - ZONE_MARGIN_WIDTH
ZONE_PREMARG_TOP <- ZONE_MARGIN_TOP - ZONE_PREMARG_WIDTH

cat(sprintf("  Margin (static): %.1f deg\n", MARGIN_THETA))

spots <- spots %>% mutate(depth_epoch = floor(FRAME / DEPTH_EPOCH_SIZE) * DEPTH_EPOCH_SIZE)

# Same-latitude depth reference for ingression cluster
ingr_lat_lo <- INGRESSION_CENTER_THETA - 10
ingr_lat_hi <- INGRESSION_CENTER_THETA + 10
same_lat_depth_ref <- spots %>%
  filter(THETA_DEG >= ingr_lat_lo & THETA_DEG <= ingr_lat_hi &
         sqrt((THETA_DEG - INGRESSION_CENTER_THETA)^2 +
              (PHI_DEG - INGRESSION_CENTER_PHI)^2) >= INGRESSION_RADIUS) %>%
  group_by(depth_epoch) %>%
  summarise(local_depth_floor = median(SPHERICAL_DEPTH, na.rm = TRUE), .groups = "drop")

spots <- spots %>% left_join(same_lat_depth_ref, by = "depth_epoch")

spots <- spots %>%
  mutate(
    dist_to_ingr = sqrt((THETA_DEG - INGRESSION_CENTER_THETA)^2 +
                        (PHI_DEG - INGRESSION_CENTER_PHI)^2),
    in_ingr_cluster = dist_to_ingr < INGRESSION_RADIUS & SPHERICAL_DEPTH >= local_depth_floor,
    zone = factor(case_when(
      THETA_DEG > MARGIN_THETA    ~ "Sub-marginal",
      THETA_DEG > ZONE_MARGIN_TOP ~ "Margin layer",
      THETA_DEG > ZONE_PREMARG_TOP ~ "Pre-marginal",
      TRUE ~ "Animal cap"
    ), levels = ZONE_LEVELS)
  )

spots_vel <- spots_vel %>%
  mutate(depth_epoch = floor(FRAME / DEPTH_EPOCH_SIZE) * DEPTH_EPOCH_SIZE) %>%
  left_join(same_lat_depth_ref, by = "depth_epoch") %>%
  mutate(
    dist_to_ingr = sqrt((THETA_DEG - INGRESSION_CENTER_THETA)^2 +
                        (PHI_DEG - INGRESSION_CENTER_PHI)^2),
    in_ingr_cluster = dist_to_ingr < INGRESSION_RADIUS & SPHERICAL_DEPTH >= local_depth_floor,
    zone = factor(case_when(
      THETA_DEG > MARGIN_THETA    ~ "Sub-marginal",
      THETA_DEG > ZONE_MARGIN_TOP ~ "Margin layer",
      THETA_DEG > ZONE_PREMARG_TOP ~ "Pre-marginal",
      TRUE ~ "Animal cap"
    ), levels = ZONE_LEVELS),
    convergence_toward_dorsal = -v_convergence * sign(PHI_DEG - INGRESSION_CENTER_PHI)
  )

# Step classification (identical algorithm to medaka)
spots_vel <- spots_vel %>%
  mutate(step_type = case_when(
    abs(v_convergence) > abs(v_epiboly) &
      convergence_toward_dorsal > 0                        ~ "Convergence",
    v_epiboly > 0                                          ~ "Epiboly",
    TRUE                                                   ~ "Animalward"
  ))

zone_summary <- spots %>% group_by(zone) %>%
  summarise(n = n(), pct = round(100 * n() / nrow(spots), 1), .groups = "drop")
cat("  Zone distribution:\n")
for (i in seq_len(nrow(zone_summary))) {
  cat(sprintf("    %-15s: %s (%.1f%%)\n", zone_summary$zone[i],
              format(zone_summary$n[i], big.mark = ","), zone_summary$pct[i]))
}

# =============================================================================
# STEP 5: PER-TRACK METRICS
# =============================================================================

cat("\n=== STEP 5: PER-TRACK METRICS ===\n")

track_metrics <- spots %>%
  arrange(TRACK_ID, FRAME) %>%
  group_by(TRACK_ID) %>%
  filter(n() >= MIN_TRACK_LENGTH) %>%
  summarise(
    n_frames    = n(),
    duration_min = (max(FRAME) - min(FRAME)) * FRAME_INTERVAL_MIN,
    start_frame = min(FRAME), end_frame = max(FRAME),
    mean_theta  = mean(THETA_DEG, na.rm = TRUE),
    mean_phi    = mean(PHI_DEG, na.rm = TRUE),
    mean_depth  = mean(SPHERICAL_DEPTH, na.rm = TRUE),
    net_dx = last(POSITION_X) - first(POSITION_X),
    net_dy = last(POSITION_Y) - first(POSITION_Y),
    net_dz = last(POSITION_Z) - first(POSITION_Z),
    net_disp   = sqrt(net_dx^2 + net_dy^2 + net_dz^2),
    total_path = sum(sqrt(diff(POSITION_X)^2 + diff(POSITION_Y)^2 + diff(POSITION_Z)^2), na.rm = TRUE),
    straightness  = ifelse(total_path > 0, net_disp / total_path, NA_real_),
    max_from_start = max(sqrt((POSITION_X - first(POSITION_X))^2 +
                              (POSITION_Y - first(POSITION_Y))^2 +
                              (POSITION_Z - first(POSITION_Z))^2)),
    confinement = ifelse(total_path > 0, max_from_start / total_path, NA_real_),
    mean_speed  = mean(inst_speed, na.rm = TRUE),
    mean_v_epiboly     = mean(v_epiboly, na.rm = TRUE),
    mean_v_convergence = mean(v_convergence, na.rm = TRUE),
    mean_conv_dorsal   = mean(-v_convergence * sign(PHI_DEG - INGRESSION_CENTER_PHI), na.rm = TRUE),
    mean_v_depth = mean(v_depth, na.rm = TRUE),
    net_dTheta = last(THETA_DEG) - first(THETA_DEG),
    net_dPhi   = last(PHI_DEG)   - first(PHI_DEG),
    mean_turning = mean(turning_angle, na.rm = TRUE),
    mean_acceleration = mean(abs(acceleration), na.rm = TRUE),
    depth_change = last(SPHERICAL_DEPTH) - first(SPHERICAL_DEPTH),
    .groups = "drop"
  ) %>%
  mutate(persistence = straightness)

track_metrics <- track_metrics %>%
  left_join(track_depth_post %>% dplyr::select(TRACK_ID, max_depth_post, n_frames_post, depth_change_post),
            by = "TRACK_ID") %>%
  mutate(
    dist_to_ingression = sqrt((mean_theta - INGRESSION_CENTER_THETA)^2 +
                              (mean_phi   - INGRESSION_CENTER_PHI)^2),
    total_v = abs(mean_v_epiboly) + abs(mean_v_convergence) + 1e-6,
    epiboly_frac     = abs(mean_v_epiboly) / total_v,
    convergence_frac = abs(mean_v_convergence) / total_v,
    net_dPhi_dorsal  = -net_dPhi * sign(mean_phi - INGRESSION_CENTER_PHI)
  )

# Kinematic ingression threshold: P90 of depth change post-onset
# Combined with spatial constraint: must also be within INGRESSION_RADIUS of cluster
post_onset_depth <- track_metrics %>%
  filter(end_frame >= INGRESSION_ONSET_FRAME & !is.na(depth_change_post))
INGRESSION_DEPTH_CHANGE_THRESH <- quantile(post_onset_depth$depth_change_post, 0.90, na.rm = TRUE)
if (is.na(INGRESSION_DEPTH_CHANGE_THRESH) || INGRESSION_DEPTH_CHANGE_THRESH <= 0) {
  INGRESSION_DEPTH_CHANGE_THRESH <- 5
}
cat(sprintf("  Ingression depth-change threshold: %.1f um (P90 of post-onset deepening)\n",
            INGRESSION_DEPTH_CHANGE_THRESH))

track_metrics <- track_metrics %>%
  mutate(
    movement_type = case_when(
      end_frame >= INGRESSION_ONSET_FRAME &
        !is.na(depth_change_post) &
        depth_change_post > INGRESSION_DEPTH_CHANGE_THRESH &
        dist_to_ingression < INGRESSION_RADIUS                  ~ "Ingression",
      abs(mean_conv_dorsal) > abs(mean_v_epiboly) &
        mean_conv_dorsal > 0                                    ~ "Convergence",
      mean_v_epiboly > 0                                        ~ "Epiboly",
      TRUE                                                      ~ "Animalward"
    ),
    zone = factor(case_when(
      mean_theta > MARGIN_THETA    ~ "Sub-marginal",
      mean_theta > ZONE_MARGIN_TOP ~ "Margin layer",
      mean_theta > ZONE_PREMARG_TOP ~ "Pre-marginal",
      TRUE ~ "Animal cap"
    ), levels = ZONE_LEVELS)
  )

# Report: do kinematic ingressors cluster near expected location?
n_ingr <- sum(track_metrics$movement_type == "Ingression")
if (n_ingr > 0) {
  ingr_tracks <- track_metrics %>% filter(movement_type == "Ingression")
  ingr_near_cluster <- sum(ingr_tracks$dist_to_ingression < INGRESSION_RADIUS)
  cat(sprintf("  Kinematic ingressors: %d total, %d (%.0f%%) near detected cluster\n",
              n_ingr, ingr_near_cluster, 100 * ingr_near_cluster / n_ingr))
  cat(sprintf("  Ingressors spatial center: theta=%.1f, phi=%.1f (cluster: theta=%.1f, phi=%.1f)\n",
              mean(ingr_tracks$mean_theta), mean(ingr_tracks$mean_phi),
              INGRESSION_CENTER_THETA, INGRESSION_CENTER_PHI))
}

cat(sprintf("  Tracks: %s (>= %d frames)\n", format(nrow(track_metrics), big.mark = ","), MIN_TRACK_LENGTH))
cat(sprintf("  Movement: Epiboly=%d, Animalward=%d, Convergence=%d, Ingression=%d\n",
            sum(track_metrics$movement_type == "Epiboly"),
            sum(track_metrics$movement_type == "Animalward"),
            sum(track_metrics$movement_type == "Convergence"),
            sum(track_metrics$movement_type == "Ingression")))

# =============================================================================
# STEP 6: MSD
# =============================================================================

cat("\n=== STEP 6: MSD ===\n")

compute_msd <- function(dt_spots, max_lag = MSD_MAX_LAG, n_sample = MSD_N_TRACKS) {
  set.seed(42)
  track_ids <- unique(dt_spots$TRACK_ID)
  if (length(track_ids) > n_sample) track_ids <- sample(track_ids, n_sample)
  dt <- as.data.table(dt_spots)[TRACK_ID %in% track_ids,
                                 .(TRACK_ID, FRAME, POSITION_X, POSITION_Y, POSITION_Z)]
  setkey(dt, TRACK_ID, FRAME)
  msd_list <- vector("list", max_lag)
  for (lag in seq_len(max_lag)) {
    d2 <- dt[, {
      n <- .N
      if (n > lag) {
        idx1 <- 1:(n - lag); idx2 <- (1 + lag):n
        valid <- (FRAME[idx2] - FRAME[idx1]) == lag
        if (sum(valid) > 0) {
          dsq <- (POSITION_X[idx2[valid]] - POSITION_X[idx1[valid]])^2 +
                 (POSITION_Y[idx2[valid]] - POSITION_Y[idx1[valid]])^2 +
                 (POSITION_Z[idx2[valid]] - POSITION_Z[idx1[valid]])^2
          list(msd = mean(dsq), n = length(dsq))
        } else list(msd = numeric(0), n = integer(0))
      } else list(msd = numeric(0), n = integer(0))
    }, by = TRACK_ID]
    if (nrow(d2) > 0 && sum(d2$n) > 0) {
      msd_list[[lag]] <- data.frame(
        lag_frames = lag, lag_min = lag * FRAME_INTERVAL_MIN,
        mean_msd = weighted.mean(d2$msd, d2$n),
        n_pairs = sum(d2$n), n_tracks = nrow(d2))
    }
  }
  bind_rows(msd_list)
}

fit_alpha <- function(msd_df) {
  d <- msd_df %>% filter(lag_frames <= 10 & mean_msd > 0)
  if (nrow(d) < 3) return(NA_real_)
  coef(lm(log10(mean_msd) ~ log10(lag_min), data = d))[2]
}

msd_global <- compute_msd(spots)
msd_global$group <- "Global"
alpha_global <- fit_alpha(msd_global)

msd_by_zone <- spots %>%
  filter(!is.na(zone)) %>%
  split(.$zone) %>%
  lapply(function(df) {
    res <- compute_msd(df, n_sample = 3000)
    res$group <- as.character(unique(df$zone)[1])
    res
  }) %>% bind_rows()

# MSD alpha over time
time_windows <- seq(0, max(spots$FRAME) - 20, by = 20)
alpha_time <- lapply(time_windows, function(t0) {
  sub <- spots %>% filter(FRAME >= t0 & FRAME <= t0 + 20)
  if (length(unique(sub$TRACK_ID)) < 100) return(NULL)
  m <- compute_msd(sub, max_lag = 10, n_sample = 2000)
  if (nrow(m) < 3) return(NULL)
  data.frame(time_min = (t0 + 10) * FRAME_INTERVAL_MIN, alpha = fit_alpha(m))
}) %>% bind_rows()

cat(sprintf("  Global alpha = %.2f\n", alpha_global))

# =============================================================================
# STEP 7: NEIGHBOR COORDINATION
# =============================================================================

cat("\n=== STEP 7: NEIGHBOR COORDINATION ===\n")

COORD_MAX_DIST <- 300  # um — capture the full correlation decay
COORD_DIST_BIN <- 5    # um — fine-grained bins for correlation function
LOCAL_RADIUS   <- 50   # um — radius for local mean velocity subtraction

set.seed(42)
sample_frames_coord <- sort(sample(unique(spots_vel$FRAME),
                                   min(30, length(unique(spots_vel$FRAME)))))

coord_results <- lapply(sample_frames_coord, function(fr) {
  sub <- spots_vel %>%
    filter(FRAME == fr, !is.na(dx), !is.na(dy), !is.na(dz)) %>%
    dplyr::select(TRACK_ID, POSITION_X, POSITION_Y, POSITION_Z, dx, dy, dz, zone)
  if (nrow(sub) < 20) return(NULL)
  if (nrow(sub) > 3000) sub <- sub %>% sample_n(3000)

  pos <- as.matrix(sub[, c("POSITION_X", "POSITION_Y", "POSITION_Z")])
  vel <- as.matrix(sub[, c("dx", "dy", "dz")])
  vel_mag <- sqrt(rowSums(vel^2))
  vel_norm <- vel / (vel_mag + 1e-10)

  # Compute local mean velocity (within LOCAL_RADIUS) for each cell
  # This separates collective flow from local neighbor coordination
  tree <- RANN::nn2(pos, pos, k = min(nrow(pos), 50), searchtype = "radius",
                    radius = LOCAL_RADIUS)
  local_vel <- t(sapply(seq_len(nrow(pos)), function(i) {
    nn_idx <- tree$nn.idx[i, ]
    nn_idx <- nn_idx[nn_idx > 0]
    colMeans(vel[nn_idx, , drop = FALSE])
  }))
  # Fluctuation velocity = individual - local mean
  vel_fluct <- vel - local_vel
  vf_mag <- sqrt(rowSums(vel_fluct^2))
  vf_norm <- vel_fluct / (vf_mag + 1e-10)

  # Shuffled baseline: permute velocity labels across cells
  shuf_idx <- sample(nrow(vel))
  vel_shuf_norm <- vel_norm[shuf_idx, ]

  n <- nrow(sub)
  n_pairs <- min(80000, n * (n - 1) / 2)
  idx_i <- sample(n, n_pairs, replace = TRUE)
  idx_j <- sample(n, n_pairs, replace = TRUE)
  valid <- idx_i != idx_j
  idx_i <- idx_i[valid]; idx_j <- idx_j[valid]
  dist_ij <- sqrt(rowSums((pos[idx_i, ] - pos[idx_j, ])^2))
  keep <- dist_ij <= COORD_MAX_DIST
  idx_i <- idx_i[keep]; idx_j <- idx_j[keep]; dist_ij <- dist_ij[keep]

  cos_raw   <- rowSums(vel_norm[idx_i, ] * vel_norm[idx_j, ])
  cos_local <- rowSums(vf_norm[idx_i, ] * vf_norm[idx_j, ])
  cos_shuf  <- rowSums(vel_shuf_norm[idx_i, ] * vel_shuf_norm[idx_j, ])

  zone_i <- as.character(sub$zone[idx_i])
  zone_j <- as.character(sub$zone[idx_j])

  data.frame(dist = dist_ij, cos_raw = cos_raw, cos_local = cos_local,
             cos_shuf = cos_shuf, zone_i = zone_i, zone_j = zone_j, frame = fr)
}) %>% bind_rows()

cat(sprintf("  %s neighbor pairs\n", format(nrow(coord_results), big.mark = ",")))

# Binned correlation functions
coord_binned <- coord_results %>%
  mutate(dist_bin = floor(dist / COORD_DIST_BIN) * COORD_DIST_BIN + COORD_DIST_BIN / 2) %>%
  group_by(dist_bin) %>%
  summarise(
    mean_cos_raw   = mean(cos_raw, na.rm = TRUE),
    se_cos_raw     = sd(cos_raw, na.rm = TRUE) / sqrt(n()),
    mean_cos_local = mean(cos_local, na.rm = TRUE),
    se_cos_local   = sd(cos_local, na.rm = TRUE) / sqrt(n()),
    mean_cos_shuf  = mean(cos_shuf, na.rm = TRUE),
    n_pairs = n(), .groups = "drop") %>%
  filter(n_pairs >= 50)

# Near-neighbor alignment by zone over time
coord_zone_time <- coord_results %>%
  filter(dist < 30, zone_i == zone_j) %>%
  mutate(time_min = frame_to_min(frame), zone = zone_i) %>%
  group_by(frame, time_min, zone) %>%
  summarise(mean_cos_raw = mean(cos_raw, na.rm = TRUE),
            mean_cos_local = mean(cos_local, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= 10)

# Distribution data for violin plots at distance ranges
coord_dist_groups <- coord_results %>%
  mutate(dist_range = cut(dist, breaks = c(0, 15, 30, 60, 100, 200, 300),
                          labels = c("<15", "15-30", "30-60", "60-100", "100-200", "200-300"))) %>%
  filter(!is.na(dist_range))

# =============================================================================
# STEP 8: THICKNESS ANALYSIS
# =============================================================================

cat("\n=== STEP 8: THICKNESS ===\n")

THICK_EPOCH <- 10
THICK_DTHETA <- 4
THICK_DPHI   <- 6

thickness_data <- spots %>%
  mutate(thick_epoch = floor(FRAME / THICK_EPOCH) * THICK_EPOCH,
         time_min_epoch = thick_epoch * FRAME_INTERVAL_MIN,
         theta_bin = floor(THETA_DEG / THICK_DTHETA) * THICK_DTHETA + THICK_DTHETA / 2,
         phi_bin   = floor(PHI_DEG / THICK_DPHI) * THICK_DPHI + THICK_DPHI / 2) %>%
  group_by(thick_epoch, time_min_epoch, theta_bin, phi_bin) %>%
  summarise(
    thickness_iqr   = IQR(SPHERICAL_DEPTH, na.rm = TRUE),
    thickness_p90   = diff(quantile(SPHERICAL_DEPTH, c(0.05, 0.95), na.rm = TRUE)),
    thickness_range = diff(range(SPHERICAL_DEPTH, na.rm = TRUE)),
    depth_mean  = mean(SPHERICAL_DEPTH, na.rm = TRUE),
    depth_sd    = sd(SPHERICAL_DEPTH, na.rm = TRUE),
    n_cells     = n(),
    .groups = "drop"
  ) %>% filter(n_cells >= 5)

thickness_zone_time <- spots %>%
  mutate(time_bin = floor(time_min / 10) * 10) %>%
  filter(!is.na(zone)) %>%
  group_by(time_bin, zone) %>%
  summarise(
    thickness_iqr = IQR(SPHERICAL_DEPTH, na.rm = TRUE),
    thickness_p90 = diff(quantile(SPHERICAL_DEPTH, c(0.05, 0.95), na.rm = TRUE)),
    thickness_range = diff(range(SPHERICAL_DEPTH, na.rm = TRUE)),
    median_depth = median(SPHERICAL_DEPTH, na.rm = TRUE),
    n_cells = n(), .groups = "drop"
  ) %>% filter(n_cells >= 20)

cat(sprintf("  Thickness bins: %d\n", nrow(thickness_data)))

# =============================================================================
# STEP 9: LOAD NUCLEAR STATS
# =============================================================================

cat("\n=== STEP 9: NUCLEAR STATS ===\n")

HAS_NUCLEAR <- dir.exists(NUCLEAR_DIR)
if (HAS_NUCLEAR) {
  nuc_time  <- read.delim(file.path(NUCLEAR_DIR, "nuclear_stats_per_time.tsv"))
  nuc_slice <- read.delim(file.path(NUCLEAR_DIR, "nuclear_stats_per_slice.tsv"))
  nuc_global <- read.delim(file.path(NUCLEAR_DIR, "nuclear_stats_global.tsv"))

  # Fix timepoints if constant (generation bug: all rows have same value)
  if (length(unique(nuc_time$timepoint)) == 1 && nrow(nuc_time) > 1) {
    cat("  Fixing per-time timepoints (all identical, assigning 0-based sequence)\n")
    nuc_time$timepoint <- seq(0, nrow(nuc_time) - 1)
  }
  if (length(unique(nuc_slice$timepoint)) == 1 && nrow(nuc_slice) > 1) {
    cat("  Fixing per-slice timepoints\n")
    n_slices <- length(unique(nuc_slice$z_slice))
    n_tp <- nrow(nuc_slice) / n_slices
    # Data is sorted by z_slice (n_tp rows per z_slice), so timepoint cycles within each z_slice block
    nuc_slice$timepoint <- rep(seq(0, n_tp - 1), times = n_slices)
  }
  cat(sprintf("  Loaded nuclear stats: %d timepoints, %d slices\n",
              nrow(nuc_time), nrow(nuc_slice)))
} else {
  cat("  No nuclear stats directory found\n")
}

# ###########################################################################
# PHASE 2: PLOTTING
# ###########################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat("  GENERATING ZEBRAFISH FIGURES\n")
cat(strrep("=", 70), "\n")

# Speed summary — all steps and moving-only (excluding zero-displacement)
speed_summary <- spots_vel %>%
  summarise(
    mean_speed = mean(inst_speed, na.rm = TRUE),
    median_speed = median(inst_speed, na.rm = TRUE),
    q95_speed = quantile(inst_speed, 0.95, na.rm = TRUE)
  )

spots_vel_moving <- spots_vel %>% filter(disp_3d > 0)
speed_summary_moving <- spots_vel_moving %>%
  summarise(
    mean_speed = mean(inst_speed, na.rm = TRUE),
    median_speed = median(inst_speed, na.rm = TRUE),
    q95_speed = quantile(inst_speed, 0.95, na.rm = TRUE)
  )

n_zero_disp <- nrow(spots_vel) - nrow(spots_vel_moving)
cat(sprintf("  Speed (all): mean=%.3f, median=%.3f um/min (N=%s)\n",
            speed_summary$mean_speed, speed_summary$median_speed, format(nrow(spots_vel), big.mark=",")))
cat(sprintf("  Speed (moving): mean=%.3f, median=%.3f um/min (N=%s, excluded %d zero-disp = %.1f%%)\n",
            speed_summary_moving$mean_speed, speed_summary_moving$median_speed,
            format(nrow(spots_vel_moving), big.mark=","), n_zero_disp, 100*n_zero_disp/nrow(spots_vel)))

# ============================================================================
# FIGURE 1: VELOCITY & ACCELERATION
# ============================================================================

cat("\n=== FIGURE 1: VELOCITY & ACCELERATION ===\n")

p1a <- ggplot(spots_vel, aes(x = inst_speed)) +
  geom_histogram(bins = 100, fill = "#2166AC", alpha = 0.4, color = NA) +
  geom_histogram(data = spots_vel_moving, bins = 100, fill = "#2166AC", alpha = 0.7, color = NA) +
  geom_vline(xintercept = speed_summary_moving$median_speed, linetype = "dashed",
             color = "#B2182B", linewidth = 0.6) +
  geom_vline(xintercept = speed_summary_moving$mean_speed, linetype = "dotted",
             color = "#B2182B", linewidth = 0.6) +
  coord_cartesian(xlim = c(0, quantile(spots_vel$inst_speed, 0.99, na.rm = TRUE))) +
  labs(title = "Instantaneous speed distribution",
       subtitle = sprintf("Mean = %.2f, Median = %.2f um/min (moving only, %d zero-disp excluded = %.1f%%).",
                          speed_summary_moving$mean_speed, speed_summary_moving$median_speed, n_zero_disp, 100*n_zero_disp/nrow(spots_vel)),
       x = "Speed (um/min)", y = "Count") +
  theme_dynamics()

p1b <- ggplot(spots_vel %>% filter(!is.na(zone)),
              aes(x = zone, y = inst_speed, fill = zone)) +
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = zone_colors) +
  coord_cartesian(ylim = c(0, quantile(spots_vel$inst_speed, 0.95, na.rm = TRUE))) +
  labs(title = "Speed by zone", x = NULL, y = "Speed (um/min)") +
  theme_dynamics() + theme(legend.position = "none")

accel_data <- spots_vel %>% filter(!is.na(acceleration)) %>% mutate(abs_accel = abs(acceleration))
ACCEL_BINS <- 60
p1c <- ggplot(accel_data, aes(x = abs_accel)) +
  geom_histogram(bins = ACCEL_BINS, fill = "#E07A3A", alpha = 0.7, color = NA) +
  geom_vline(xintercept = median(accel_data$abs_accel, na.rm = TRUE),
             linetype = "dashed", color = "#B2182B", linewidth = 0.6) +
  geom_vline(xintercept = mean(accel_data$abs_accel, na.rm = TRUE),
             linetype = "dotted", color = "#2166AC", linewidth = 0.6) +
  coord_cartesian(xlim = c(0, quantile(accel_data$abs_accel, 0.99, na.rm = TRUE))) +
  labs(title = "Acceleration magnitude distribution",
       subtitle = sprintf("Median = %.3f, Mean = %.3f um/min² (n=%s, %d bins)",
                          median(accel_data$abs_accel, na.rm = TRUE),
                          mean(accel_data$abs_accel, na.rm = TRUE),
                          format(nrow(accel_data), big.mark = ","), ACCEL_BINS),
       x = "|Acceleration| (um/min²)", y = "Count") +
  theme_dynamics()

# Turning angle — smoothing sensitivity validation
smooth_colors <- c("#D73027", "#FC8D59", "#91BFDB", "#4575B4")
names(smooth_colors) <- SMOOTH_LABELS

smooth_median_df <- as.data.frame(smooth_val_summary)

p1d <- ggplot(as.data.frame(smooth_val_results),
              aes(x = turning_angle, color = smooth_level, fill = smooth_level)) +
  geom_density(alpha = 0.12, linewidth = 0.7) +
  geom_vline(data = smooth_median_df,
             aes(xintercept = median_turn, color = smooth_level),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  scale_color_manual(values = smooth_colors, name = "Smoothing") +
  scale_fill_manual(values = smooth_colors, name = "Smoothing") +
  labs(title = "Turning angle: smoothing sensitivity test",
       subtitle = sprintf("Median turning angle -- raw: %.1f deg, current (k=%d): %.1f deg, aggressive (k=%d): %.1f deg (%+.1f%%)",
                          smooth_val_summary[smooth_k == 1, median_turn],
                          SMOOTH_K,
                          smooth_val_summary[smooth_k == SMOOTH_K, median_turn],
                          max(SMOOTH_LEVELS),
                          smooth_val_summary[smooth_k == max(SMOOTH_LEVELS), median_turn],
                          aggr_turn_pct),
       x = "Turning angle (deg)", y = "Density") +
  theme_dynamics()

# Speed stability across smoothing levels (control)
p1h <- ggplot(as.data.frame(smooth_val_results),
              aes(x = inst_speed, color = smooth_level, fill = smooth_level)) +
  geom_density(alpha = 0.12, linewidth = 0.7) +
  coord_cartesian(xlim = c(0, quantile(smooth_val_results$inst_speed, 0.99, na.rm = TRUE))) +
  scale_color_manual(values = smooth_colors, name = "Smoothing") +
  scale_fill_manual(values = smooth_colors, name = "Smoothing") +
  labs(title = "Speed stability across smoothing levels",
       subtitle = sprintf("Median speed -- raw: %.3f, aggressive (k=%d): %.3f um/min (%+.1f%%)",
                          smooth_val_summary[smooth_k == 1, median_speed],
                          max(SMOOTH_LEVELS),
                          smooth_val_summary[smooth_k == max(SMOOTH_LEVELS), median_speed],
                          aggr_speed_pct),
       x = "Speed (um/min)", y = "Density") +
  theme_dynamics()

# Speed over time
speed_time <- spots_vel %>%
  filter(!is.na(zone)) %>%
  mutate(time_bin = floor(time_min / 20) * 20) %>%
  group_by(time_bin, zone) %>%
  summarise(mean_speed = mean(inst_speed, na.rm = TRUE),
            se = sd(inst_speed, na.rm = TRUE) / sqrt(n()),
            n = n(), .groups = "drop") %>%
  filter(n >= 30)

p1e <- ggplot(speed_time, aes(x = time_bin, y = mean_speed, color = zone)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = mean_speed - se, ymax = mean_speed + se, fill = zone),
              alpha = 0.12, color = NA) +
  scale_color_manual(values = zone_colors) +
  scale_fill_manual(values = zone_colors) +
  labs(title = "Mean speed over time by zone", x = "Time (min)", y = "Speed (um/min)") +
  guides(fill = "none") + theme_dynamics()

# Velocity components by zone over time
comp_time <- spots_vel %>%
  filter(!is.na(zone)) %>%
  mutate(time_bin = floor(time_min / 20) * 20) %>%
  group_by(time_bin, zone) %>%
  summarise(v_epiboly = mean(v_epiboly, na.rm = TRUE),
            v_convergence = mean(v_convergence, na.rm = TRUE),
            v_depth = mean(v_depth, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= 30) %>%
  pivot_longer(cols = c(v_epiboly, v_convergence, v_depth),
               names_to = "component", values_to = "velocity")

p1f <- ggplot(comp_time, aes(x = time_bin, y = velocity, color = component)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  facet_wrap(~zone, scales = "free_y") +
  scale_color_manual(values = c("v_epiboly" = "#EF8A62", "v_convergence" = "#1B7837",
                                "v_depth" = "#9970AB"),
                     labels = c("Epiboly", "Convergence", "Depth"), name = "Component") +
  labs(title = "Velocity components over time by zone",
       x = "Time (min)", y = "Velocity (um/min)") +
  theme_dynamics(8)

# Velocity components by latitude
speed_theta <- spots_vel %>%
  mutate(theta_bin = floor(THETA_DEG / THETA_BIN) * THETA_BIN + THETA_BIN / 2) %>%
  group_by(theta_bin) %>%
  summarise(speed = mean(inst_speed, na.rm = TRUE),
            v_epiboly = mean(v_epiboly, na.rm = TRUE),
            v_convergence = mean(v_convergence, na.rm = TRUE),
            v_depth = mean(v_depth, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= FLOW_MIN_N) %>%
  pivot_longer(cols = c(speed, v_epiboly, v_convergence, v_depth),
               names_to = "component", values_to = "velocity")

p1g <- ggplot(speed_theta, aes(x = theta_bin, y = velocity, color = component)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  scale_color_manual(values = c("speed" = "black", "v_epiboly" = "#EF8A62",
                                "v_convergence" = "#1B7837", "v_depth" = "#9970AB"),
                     name = "Component") +
  labs(title = "Velocity components by latitude",
       x = THETA_LABEL, y = "Velocity (um/min)") +
  theme_dynamics()

fig1 <- (p1a + p1b) / (p1c + p1d) / (p1h + p1e) / (p1f + p1g) +
  plot_layout(heights = c(1, 1, 1, 1)) +
  plot_annotation(
    title = sprintf("%s -- Velocity, Acceleration & Smoothing Validation", SPECIES),
    theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14))
  )
save_pdf(fig1, "zebrafish_dynamics_01_velocity.pdf", width = 16, height = 22)

# ============================================================================
# FIGURE 2: TRACK STRAIGHTNESS & PERSISTENCE
# ============================================================================

cat("\n=== FIGURE 2: TRACK PROPERTIES ===\n")

p2a <- ggplot(track_metrics %>% filter(!is.na(zone)),
              aes(x = zone, y = straightness, fill = zone)) +
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = zone_colors) +
  labs(title = "Track straightness by zone", x = NULL, y = "Straightness") +
  theme_dynamics() + theme(legend.position = "none")

p2b <- ggplot(track_metrics %>% filter(!is.na(mean_turning), !is.na(zone)),
              aes(x = zone, y = mean_turning, fill = zone)) +
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = zone_colors) +
  labs(title = "Mean turning angle by zone", x = NULL, y = "Turning angle (deg)") +
  theme_dynamics() + theme(legend.position = "none")

p2c <- ggplot(track_metrics, aes(x = movement_type, y = confinement, fill = movement_type)) +
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = mvt_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Confinement by movement type", x = NULL, y = "Confinement ratio",
       subtitle = "Max displacement / total path. 1 = straight line, 0 = returns to start.") +
  theme_dynamics() + theme(legend.position = "none")

p2d <- ggplot(track_metrics, aes(x = mean_speed, y = straightness, color = movement_type)) +
  geom_point(alpha = 0.15, size = 0.5) +
  scale_color_manual(values = mvt_colors) +
  coord_cartesian(xlim = c(0, quantile(track_metrics$mean_speed, 0.99, na.rm = TRUE))) +
  labs(title = "Straightness vs speed", x = "Mean speed (um/min)", y = "Straightness") +
  theme_dynamics()

fig2 <- (p2a + p2b) / (p2c + p2d) +
  plot_annotation(
    title = sprintf("%s — Track Straightness, Persistence & Confinement", SPECIES),
    subtitle = sprintf("N = %s tracks", format(nrow(track_metrics), big.mark = ",")),
    theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14),
                  plot.subtitle = element_text(hjust = 0, size = 10, color = "grey50"))
  )
save_pdf(fig2, "zebrafish_dynamics_02_track_properties.pdf", width = 14, height = 10)

# ============================================================================
# FIGURE 3: NEIGHBOR COORDINATION
# ============================================================================

cat("\n=== FIGURE 3: NEIGHBOR COORDINATION ===\n")

# 3a: Velocity correlation function — raw + local + shuffled baseline
corr_long <- coord_binned %>%
  dplyr::select(dist_bin, mean_cos_raw, mean_cos_local, mean_cos_shuf) %>%
  pivot_longer(-dist_bin, names_to = "type", values_to = "cos") %>%
  mutate(type = recode(type,
    "mean_cos_raw"   = "Total (raw)",
    "mean_cos_local" = "Local (collective flow subtracted)",
    "mean_cos_shuf"  = "Shuffled (random expectation)"))

# Find correlation length (where local curve crosses 0)
corr_len <- coord_binned %>%
  filter(mean_cos_local > 0) %>%
  summarise(corr_length = max(dist_bin)) %>%
  pull(corr_length)

p3a <- ggplot(corr_long, aes(x = dist_bin, y = cos, color = type, linewidth = type)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = corr_len, linetype = "dotted", color = "#E08214", linewidth = 0.5) +
  annotate("text", x = corr_len + 5, y = max(coord_binned$mean_cos_raw) * 0.9,
           label = sprintf("Corr. length ~%.0f um", corr_len),
           hjust = 0, size = 3, color = "#E08214") +
  scale_color_manual(values = c("Total (raw)" = "#2166AC",
                                "Local (collective flow subtracted)" = "#D73027",
                                "Shuffled (random expectation)" = "grey60"),
                     name = NULL) +
  scale_linewidth_manual(values = c("Total (raw)" = 0.9,
                                    "Local (collective flow subtracted)" = 0.9,
                                    "Shuffled (random expectation)" = 0.5),
                         guide = "none") +
  labs(title = "Velocity correlation vs distance",
       subtitle = "1 = identical direction, 0 = uncorrelated, -1 = opposite",
       x = "Distance between cells (um)", y = "Mean cosine similarity") +
  theme_dynamics() +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# 3b: Distribution of cosine similarity at different distances (violin + box)
p3b <- ggplot(coord_dist_groups, aes(x = dist_range, y = cos_raw, fill = dist_range)) +
  geom_violin(alpha = 0.5, color = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.1, fill = "white", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  labs(title = "Alignment distribution by distance",
       subtitle = "Nearby cells: tight alignment; distant cells: broad spread",
       x = "Distance range (um)", y = "Cosine similarity") +
  theme_dynamics() + theme(legend.position = "none")

# 3c: Local coordination correlation function — zoomed to show structure
# after subtracting collective flow
p3c <- ggplot(coord_binned, aes(x = dist_bin, y = mean_cos_local)) +
  geom_ribbon(aes(ymin = mean_cos_local - se_cos_local,
                  ymax = mean_cos_local + se_cos_local),
              alpha = 0.2, fill = "#D73027") +
  geom_line(linewidth = 0.9, color = "#D73027") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_cartesian(xlim = c(0, min(200, max(coord_binned$dist_bin)))) +
  labs(title = "Local coordination (collective flow removed)",
       subtitle = "True neighbor-neighbor coordination beyond shared drift",
       x = "Distance between cells (um)", y = "Mean cosine similarity") +
  theme_dynamics()

# 3d: Near-neighbor alignment by zone over time
if (nrow(coord_zone_time) > 0 && length(unique(coord_zone_time$zone)) > 1) {
  p3d <- ggplot(coord_zone_time %>% filter(!is.na(zone)),
                aes(x = time_min, y = mean_cos_raw, color = zone)) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = zone_colors, name = "Zone") +
    labs(title = "Near-neighbor alignment (<30 um) by zone",
         subtitle = "How coordination varies across different embryo regions",
         x = "Time (min)", y = "Mean cosine similarity") +
    theme_dynamics()
} else {
  p3d <- ggplot(coord_binned, aes(x = dist_bin, y = n_pairs)) +
    geom_col(fill = "#2166AC", alpha = 0.5) +
    labs(title = "Number of pairs per distance bin", x = "Distance (um)", y = "N pairs") +
    theme_dynamics()
}

fig3 <- (p3a + p3b) / (p3c + p3d) +
  plot_annotation(
    title = sprintf("%s — Neighbor Velocity Coordination", SPECIES),
    subtitle = "How aligned are velocity directions between nearby cells?",
    theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14),
                  plot.subtitle = element_text(hjust = 0, size = 10, color = "grey30"))
  )
save_pdf(fig3, "zebrafish_dynamics_03_neighbors.pdf", width = 14, height = 10)

# ============================================================================
# FIGURE 6: DIRECTIONAL DECOMPOSITION
# ============================================================================

cat("\n=== FIGURE 6: DIRECTIONAL DECOMPOSITION ===\n")

dir_angle <- spots_vel %>%
  mutate(angle = atan2(v_convergence, v_epiboly) * 180 / pi)

p4a <- ggplot(dir_angle, aes(x = angle)) +
  geom_histogram(bins = 72, fill = "#2166AC", alpha = 0.7, color = NA) +
  scale_x_continuous(breaks = c(-180, -90, 0, 90, 180),
                     labels = c("Div.", "Anti-epi.", "Epi.", "Conv.", "Div.")) +
  labs(title = "Direction of tangential movement", x = "Direction (deg)", y = "Count") +
  theme_dynamics()

p4b <- ggplot(spots_vel, aes(x = step_type, y = inst_speed, fill = step_type)) +
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = step_colors) +
  coord_cartesian(ylim = c(0, quantile(spots_vel$inst_speed, 0.95, na.rm = TRUE))) +
  labs(title = "Speed by direction", x = NULL, y = "Speed (um/min)") +
  theme_dynamics() + theme(axis.text.x = element_text(angle = 25, hjust = 1),
                           legend.position = "none")

step_time <- spots_vel %>%
  mutate(time_bin = floor(time_min / 20) * 20) %>%
  group_by(time_bin, step_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(time_bin) %>%
  mutate(pct = 100 * n / sum(n))

p4c <- ggplot(step_time, aes(x = time_bin, y = pct, fill = step_type)) +
  geom_area(alpha = 0.8, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = step_colors) +
  labs(title = "Movement direction over time", x = "Time (min)", y = "%", fill = "Direction") +
  theme_dynamics()

mvt_prop <- track_metrics %>% filter(!is.na(zone)) %>%
  group_by(zone, movement_type) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(zone) %>% mutate(pct = 100 * n / sum(n))

p4d <- ggplot(mvt_prop, aes(x = zone, y = pct, fill = movement_type)) +
  geom_col(position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = mvt_colors) +
  labs(title = "Track movement type by zone", x = NULL, y = "%", fill = "Type") +
  theme_dynamics()

fig4 <- (p4a + p4b) / (p4c + p4d) +
  plot_annotation(
    title = sprintf("%s — Directional Decomposition", SPECIES),
    theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14))
  )
save_pdf(fig4, "zebrafish_dynamics_06_direction.pdf", width = 14, height = 10)

# ============================================================================
# FIGURE 7: MSD
# ============================================================================

cat("\n=== FIGURE 7: MSD ===\n")

p5a <- ggplot(msd_global, aes(x = lag_min, y = mean_msd)) +
  geom_line(linewidth = 1, color = "#2166AC") +
  geom_point(size = 1.5, color = "#2166AC") +
  geom_abline(slope = 1, intercept = log10(msd_global$mean_msd[1]) - log10(msd_global$lag_min[1]),
              linetype = "dashed", color = "grey50") +
  geom_abline(slope = 2, intercept = log10(msd_global$mean_msd[1]) - 2*log10(msd_global$lag_min[1]),
              linetype = "dotted", color = "grey50") +
  scale_x_log10() + scale_y_log10() + annotation_logticks(sides = "bl") +
  labs(title = "Global MSD",
       subtitle = sprintf("alpha = %.2f", alpha_global),
       x = "Time lag (min)", y = expression(MSD~(mu*m^2))) +
  theme_dynamics()

p5b <- ggplot(msd_by_zone, aes(x = lag_min, y = mean_msd, color = group)) +
  geom_line(linewidth = 0.8) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = zone_colors) +
  annotation_logticks(sides = "bl") +
  labs(title = "MSD by zone", x = "Time lag (min)", y = expression(MSD~(mu*m^2)), color = "Zone") +
  theme_dynamics()

p5c <- ggplot(alpha_time %>% filter(!is.na(alpha)), aes(x = time_min, y = alpha)) +
  geom_line(linewidth = 0.8, color = "#B2182B") +
  geom_point(size = 1.5, color = "#B2182B") +
  geom_hline(yintercept = c(1, 2), linetype = c("dashed", "dotted"), color = "grey50") +
  labs(title = "MSD exponent over time", x = "Time (min)", y = "alpha") +
  theme_dynamics()

p5d <- ggplot(msd_global %>% filter(lag_frames <= 10),
              aes(x = lag_min, y = mean_msd)) +
  geom_line(linewidth = 0.8, color = "#2166AC") +
  geom_point(size = 1.5, color = "#2166AC") +
  labs(title = "MSD short-lag detail (first 10 lags)",
       x = "Time lag (min)", y = expression(MSD~(mu*m^2))) +
  theme_dynamics()

fig5 <- (p5a + p5b) / (p5c + p5d) +
  plot_annotation(
    title = sprintf("%s — Mean Square Displacement", SPECIES),
    theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14))
  )
save_pdf(fig5, "zebrafish_dynamics_07_msd.pdf", width = 14, height = 10)

# ============================================================================
# FIGURE 8: THICKNESS
# ============================================================================

cat("\n=== FIGURE 8: THICKNESS ===\n")

p6a <- ggplot(thickness_zone_time, aes(x = time_bin, y = thickness_p90, color = zone)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = zone_colors) +
  labs(title = "Blastoderm thickness over time (P5-P95 depth range)",
       x = "Time (min)", y = "Thickness P5-P95 (um)", color = "Zone") +
  theme_dynamics()

p6b <- ggplot(spots %>% filter(!is.na(zone)) %>%
                mutate(time_bin = floor(time_min / 20) * 20) %>%
                group_by(time_bin, zone) %>%
                summarise(median_depth = median(SPHERICAL_DEPTH, na.rm = TRUE),
                          n = n(), .groups = "drop") %>%
                filter(n >= 20),
              aes(x = time_bin, y = median_depth, color = zone)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = zone_colors) +
  labs(title = "Median depth over time by zone",
       x = "Time (min)", y = "Depth (um)", color = "Zone") +
  theme_dynamics()

# Depth distribution by zone
p6c <- ggplot(spots %>% filter(!is.na(zone)),
              aes(x = zone, y = SPHERICAL_DEPTH, fill = zone)) +
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = zone_colors) +
  labs(title = "Depth distribution by zone", x = NULL, y = "Depth (um)") +
  theme_dynamics() + theme(legend.position = "none")

# Estimated cell layers
# Cell spacing = NN distance (accounts for nuclear diameter + inter-nuclear gap)
if (HAS_NUCLEAR) {
  nuc_spacing <- nuc_time %>%
    mutate(time_min = timepoint * FRAME_INTERVAL_MIN,
           time_bin = floor(time_min / 10) * 10,
           cell_spacing = nn3d_mean_um) %>%
    group_by(time_bin) %>%
    summarise(cell_spacing = mean(cell_spacing, na.rm = TRUE),
              nuc_diam = mean(eqsph_diam_mean_um, na.rm = TRUE),
              nn_dist  = mean(nn3d_mean_um, na.rm = TRUE), .groups = "drop")
  cell_layers <- thickness_zone_time %>%
    left_join(nuc_spacing, by = "time_bin") %>%
    mutate(cell_spacing = ifelse(is.na(cell_spacing), mean(nuc_spacing$cell_spacing), cell_spacing),
           n_layers = thickness_p90 / cell_spacing)
  spacing_label <- sprintf("cell spacing = NN dist (mean %.1f um)", mean(nuc_spacing$cell_spacing))
} else {
  mean_cell_spacing <- 12
  cell_layers <- thickness_zone_time %>%
    mutate(n_layers = thickness_p90 / mean_cell_spacing)
  spacing_label <- sprintf("cell spacing = %.0f um (fallback)", mean_cell_spacing)
}

p6d <- ggplot(cell_layers, aes(x = time_bin, y = n_layers, color = zone)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = zone_colors) +
  labs(title = "Estimated cell layers (thickness / NN distance)",
       subtitle = spacing_label,
       x = "Time (min)", y = "Layers", color = "Zone") +
  theme_dynamics()

fig6 <- (p6a + p6b) / (p6c + p6d) +
  plot_annotation(
    title = sprintf("%s — Blastoderm Thickness & Depth", SPECIES),
    theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14))
  )
save_pdf(fig6, "zebrafish_dynamics_08_thickness.pdf", width = 14, height = 12)

# ============================================================================
# ============================================================================
# FIGURE 9: NUCLEAR SPATIAL STRUCTURE (spatial evidence first)
# Each row pairs a heatmap (time x z, left) with its z-profile (right).
# This figure establishes WHY superficial z-slices have inflated metrics.
# ============================================================================

cat("\n=== FIGURE 9: NUCLEAR SPATIAL STRUCTURE ===\n")

if (HAS_NUCLEAR) {

  # --- Detect and flag bad timepoints (segmentation failures) ---
  nuc_time$bad_tp <- FALSE
  for (iter in 1:5) {
    good_idx <- which(!nuc_time$bad_tp)
    lo <- loess(n_nuclei_3d ~ timepoint, data = nuc_time[good_idx, ], span = 0.3)
    expected <- predict(lo, newdata = nuc_time)
    nuc_time$bad_tp <- !is.na(expected) & (nuc_time$n_nuclei_3d < 0.75 * expected)
  }

  bad_tps <- nuc_time$timepoint[nuc_time$bad_tp]
  cat(sprintf("  Bad timepoints (segmentation failures): %d detected [%s]\n",
              length(bad_tps),
              if (length(bad_tps) > 0) paste(range(bad_tps), collapse = "-") else "none"))
  nuc_time_clean <- nuc_time %>% filter(!bad_tp)

  # Build bad-tp shading rectangles for annotating plots
  if (length(bad_tps) > 0) {
    bad_runs <- rle(nuc_time$bad_tp)
    run_end <- cumsum(bad_runs$lengths)
    run_start <- c(1, head(run_end, -1) + 1)
    bad_rects <- data.frame(
      xmin = nuc_time$timepoint[run_start[bad_runs$values]],
      xmax = nuc_time$timepoint[run_end[bad_runs$values]]
    )
  } else {
    bad_rects <- data.frame(xmin = numeric(0), xmax = numeric(0))
  }
  add_bad_tp_shading <- function(p) {
    if (nrow(bad_rects) > 0) {
      p <- p + geom_rect(data = bad_rects, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                         fill = "grey80", alpha = 0.3, inherit.aes = FALSE)
    }
    p
  }

  # ── Shared depth-zone infrastructure (used by Figs 09 & 10) ──────────────
  # Two groups: Edge (outside core z-range) and Core (within core z-range).
  max_z    <- max(nuc_slice$z_slice, na.rm = TRUE) + 1

  nuc_slice_binned <- nuc_slice %>%
    filter(n_nuclei > 0, !timepoint %in% bad_tps) %>%
    left_join(nuc_time %>% dplyr::select(timepoint, z_core_lo, z_core_hi),
              by = "timepoint") %>%
    mutate(z_bin = ifelse(!is.na(z_core_lo) &
                          z_slice >= z_core_lo & z_slice <= z_core_hi,
                          "Core", "Edge")) %>%
    mutate(z_bin = factor(z_bin, levels = c("Edge", "Core")))

  # 2D diameter by depth zone (n-weighted)
  z_strat_diam <- nuc_slice_binned %>%
    filter(!is.na(eqcircle_diam_mean_um)) %>%
    group_by(timepoint, z_bin) %>%
    summarise(diam_nwt = weighted.mean(eqcircle_diam_mean_um, n_nuclei, na.rm = TRUE),
              n_total  = sum(n_nuclei), .groups = "drop")

  # 2D NN distance by depth zone (n-weighted)
  z_strat_nn <- nuc_slice_binned %>%
    filter(!is.na(nn_dist_mean_um)) %>%
    group_by(timepoint, z_bin) %>%
    summarise(nn_nwt = weighted.mean(nn_dist_mean_um, n_nuclei, na.rm = TRUE),
              .groups = "drop")

  # Monocolor depth palette: light = Edge, dark = Core
  depth_levels <- c("Edge", "Core")
  zbin_pal <- c("Edge"  = "#93C5FD",  # light blue
                "Core"  = "#1E40AF",  # dark blue
                "3D equiv-sphere" = "#8B5CF6",
                "3D NN"           = "#8B5CF6",
                "3D bbox major"   = "#E07A3A",
                "3D bbox minor"   = "#F59E0B")

  # ── Core-z boundaries (slices with >25 % of peak nuclei count) ───────────
  has_core_z <- "z_core_lo" %in% names(nuc_time) && !all(is.na(nuc_time$z_core_lo))
  if (has_core_z) {
    core_bounds <- nuc_time %>% filter(!bad_tp) %>%
      dplyr::select(timepoint, z_core_lo, z_core_hi) %>% filter(!is.na(z_core_lo))
  }
  add_core_lines <- function(p) {
    if (has_core_z) {
      p <- p +
        geom_line(data = core_bounds, aes(x = timepoint, y = z_core_lo),
                  color = "white", linewidth = 0.4, linetype = "dashed", inherit.aes = FALSE) +
        geom_line(data = core_bounds, aes(x = timepoint, y = z_core_hi),
                  color = "white", linewidth = 0.4, linetype = "dashed", inherit.aes = FALSE)
    }
    p
  }

  # ── Data subsets for heatmaps and z-profiles ─────────────────────────────
  nuc_slice_clean <- nuc_slice %>% filter(!timepoint %in% bad_tps)
  nuc_slice_nz    <- nuc_slice_clean %>% filter(n_nuclei > 0)

  all_tps  <- sort(unique(nuc_slice_clean$timepoint))
  pick_idx <- round(seq(1, length(all_tps), length.out = min(8, length(all_tps))))
  pick_tps <- all_tps[pick_idx]
  df_zprof <- nuc_slice_nz %>%
    filter(timepoint %in% pick_tps) %>%
    mutate(tp_label = factor(paste0("t=", timepoint), levels = paste0("t=", pick_tps)))
  tp_pal <- scales::viridis_pal(option = "viridis")(length(pick_tps))

  # Row 1: Density heatmap + Nuclei count z-profile
  p9a <- ggplot(nuc_slice_nz, aes(timepoint, z_slice, fill = density_per_um2)) +
    geom_tile() + scale_fill_viridis_c(option = "inferno", name = expression("n/"*mu*m^2)) +
    labs(title = "Nuclear density per z-slice", x = "Timepoint", y = "Z slice",
         subtitle = "Dashed white = core z-range (>25% peak nuclei)") + theme_dynamics()
  p9a <- add_core_lines(p9a)

  p9b <- ggplot(df_zprof, aes(z_slice, n_nuclei, color = tp_label)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = setNames(tp_pal, levels(df_zprof$tp_label)), name = "Timepoint") +
    labs(title = "Nuclei count z-profile", x = "Z slice", y = "Nuclei per slice",
         subtitle = "Superficial slices have fewer nuclei; most are in the core") + theme_dynamics()

  # Row 2: Equiv-circle diameter heatmap + Diameter z-profile (KEY evidence)
  p9c <- ggplot(nuc_slice_nz, aes(timepoint, z_slice, fill = eqcircle_diam_median_um)) +
    geom_tile() + scale_fill_viridis_c(option = "cividis", name = expression(mu*m)) +
    labs(title = "Equiv-circle diameter (2D median)", x = "Timepoint", y = "Z slice") + theme_dynamics()
  p9c <- add_core_lines(p9c)

  p9d <- ggplot(df_zprof, aes(z_slice, eqcircle_diam_mean_um, color = tp_label)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = setNames(tp_pal, levels(df_zprof$tp_label)), name = "Timepoint") +
    labs(title = "Equiv-circle diameter z-profile",
         subtitle = "Edge slices show inflated & variable diameters; core may differ by species",
         x = "Z slice", y = expression("Mean diam ("*mu*m*")")) + theme_dynamics()

  # Row 3: Nuclear area heatmap + Area z-profile
  p9e <- ggplot(nuc_slice_nz, aes(timepoint, z_slice, fill = size_median_um2)) +
    geom_tile() + scale_fill_viridis_c(option = "magma", name = expression(mu*m^2)) +
    labs(title = "Nuclear area (2D median)", x = "Timepoint", y = "Z slice") + theme_dynamics()
  p9e <- add_core_lines(p9e)

  p9f <- ggplot(df_zprof, aes(z_slice, size_mean_um2, color = tp_label)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = setNames(tp_pal, levels(df_zprof$tp_label)), name = "Timepoint") +
    labs(title = "Nuclear area z-profile", x = "Z slice",
         y = expression("Mean area ("*mu*m^2*")")) + theme_dynamics()

  # Row 4: NN distance heatmap + NN distance z-profile
  p9g <- ggplot(nuc_slice_nz, aes(timepoint, z_slice, fill = nn_dist_median_um)) +
    geom_tile() + scale_fill_viridis_c(option = "mako", name = expression(mu*m)) +
    labs(title = "NN distance (2D median)", x = "Timepoint", y = "Z slice") + theme_dynamics()
  p9g <- add_core_lines(p9g)

  p9h <- ggplot(df_zprof, aes(z_slice, nn_dist_mean_um, color = tp_label)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = setNames(tp_pal, levels(df_zprof$tp_label)), name = "Timepoint") +
    labs(title = "NN distance z-profile", x = "Z slice",
         y = expression("Mean NN dist ("*mu*m*")")) + theme_dynamics()

  fig9 <- (p9a + p9b) / (p9c + p9d) / (p9e + p9f) / (p9g + p9h) +
    plot_annotation(
      title = sprintf("%s -- Nuclear Spatial Structure", SPECIES),
      subtitle = sprintf("Left: heatmaps (time x z-slice). Right: z-profiles at 8 timepoints. %d bad timepoints excluded.",
                         length(bad_tps)),
      theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14),
                    plot.subtitle = element_text(hjust = 0, size = 10, color = "grey50"))
    )
  save_pdf(fig9, "zebrafish_dynamics_09_nuclear_spatial.pdf", width = 14, height = 20)

  # ==========================================================================
  # FIGURE 10: NUCLEAR TEMPORAL TRENDS (depth-filtered)
  # 2D metrics stratified into Edge vs Core using per-timepoint core bounds.
  # 3D metrics shown for all nuclei and core-filtered (if columns exist).
  # ==========================================================================
  cat("\n=== FIGURE 10: NUCLEAR TEMPORAL TRENDS ===\n")

  nuc_time_clean <- nuc_time_clean %>%
    mutate(aspect_ratio = bbox_major_mean_um / bbox_minor_mean_um)

  # Check whether core-filtered 3D columns exist in per_time.tsv
  has_core_3d <- "core_vol_median_um3" %in% names(nuc_time)
  core_3d_label <- if (has_core_3d) " (core-filtered)" else " (all nuclei, rerun nuclear_stats.py for core filter)"

  # Row 1: Nuclear count (with bad-tp shading) + Nuclei fraction by depth zone
  p10a <- ggplot(nuc_time, aes(timepoint, n_nuclei_3d)) +
    geom_line(linewidth = 0.5, color = "grey30") +
    labs(title = "Nuclei count (3D)", x = "Timepoint", y = "N nuclei",
         subtitle = sprintf("Grey bands = %d bad timepoints (segmentation failures)", length(bad_tps))) +
    theme_dynamics()
  p10a <- add_bad_tp_shading(p10a)

  z_frac <- z_strat_diam %>%
    group_by(timepoint) %>% mutate(frac = n_total / sum(n_total)) %>% ungroup()

  p10b <- ggplot(z_frac, aes(timepoint, frac, fill = z_bin)) +
    geom_area(alpha = 0.8) +
    scale_fill_manual(values = zbin_pal[depth_levels], name = "Depth zone") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Nuclei depth distribution over time",
         subtitle = "Fraction of 2D nuclei in Edge vs Core zones",
         x = "Timepoint", y = "Fraction of nuclei") +
    theme_dynamics()

  # Row 2: 2D diameter by depth zone (+ 3D equiv-sphere) + 3D volume
  z_strat_3d_diam <- nuc_time_clean %>%
    transmute(timepoint, z_bin = "3D equiv-sphere",
              diam_nwt = if (has_core_3d) core_eqsph_diam_median_um else eqsph_diam_median_um,
              n_total = n_nuclei_3d)
  df_diam_strat <- bind_rows(z_strat_diam, z_strat_3d_diam)

  p10c <- ggplot(df_diam_strat, aes(timepoint, diam_nwt, color = z_bin)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = zbin_pal, name = "Depth zone") +
    labs(title = "Nuclear diameter by depth zone",
         subtitle = paste0("2D n-weighted diam (Edge vs Core) + 3D equiv-sphere", core_3d_label),
         x = "Timepoint", y = expression("Diameter ("*mu*m*")")) +
    theme_dynamics()

  vol_median_col <- if (has_core_3d) nuc_time_clean$core_vol_median_um3 else nuc_time_clean$vol_median_um3
  vol_mean_col   <- if (has_core_3d) nuc_time_clean$core_vol_mean_um3   else nuc_time_clean$vol_mean_um3
  p10d <- ggplot(nuc_time_clean, aes(timepoint)) +
    geom_line(aes(y = vol_median_col), linewidth = 0.5, color = "#E07A3A") +
    geom_line(aes(y = vol_mean_col), linewidth = 0.3, color = "#E07A3A",
              linetype = "dashed", alpha = 0.6) +
    labs(title = paste0("Nuclear volume (solid=median, dashed=mean)", core_3d_label),
         x = "Timepoint", y = expression("Volume ("*mu*m^3*")")) +
    theme_dynamics()

  # Row 3: 2D NN distance by depth zone (+ 3D NN) + 3D bbox major & minor axes
  nn3d_col <- if (has_core_3d) nuc_time_clean$core_nn3d_median_um else nuc_time_clean$nn3d_median_um
  z_strat_3d_nn <- nuc_time_clean %>%
    transmute(timepoint, z_bin = "3D NN", nn_nwt = nn3d_col)
  df_nn_strat <- bind_rows(z_strat_nn, z_strat_3d_nn)

  p10e <- ggplot(df_nn_strat, aes(timepoint, nn_nwt, color = z_bin)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = zbin_pal, name = "Depth zone") +
    labs(title = "Internuclear distance by depth zone",
         subtitle = paste0("2D n-weighted NN (Edge vs Core) + 3D NN", core_3d_label),
         x = "Timepoint", y = expression("NN distance ("*mu*m*")")) +
    theme_dynamics()

  # Bbox major & minor axes (3D bounding-box longest & shortest dimensions)
  bbox_major_col <- if (has_core_3d) nuc_time_clean$core_bbox_major_mean_um else nuc_time_clean$bbox_major_mean_um
  bbox_minor_col <- if (has_core_3d) nuc_time_clean$core_bbox_minor_mean_um else nuc_time_clean$bbox_minor_mean_um
  df_bbox <- tibble(
    timepoint = rep(nuc_time_clean$timepoint, 2),
    axis      = rep(c("3D bbox major", "3D bbox minor"), each = nrow(nuc_time_clean)),
    value     = c(bbox_major_col, bbox_minor_col)
  )

  p10f <- ggplot(df_bbox, aes(timepoint, value, color = axis)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = zbin_pal[c("3D bbox major", "3D bbox minor")], name = "Axis") +
    labs(title = paste0("Nuclear bbox axes (major & minor)", core_3d_label),
         x = "Timepoint", y = expression("Axis length ("*mu*m*")")) +
    theme_dynamics()

  fig10 <- (p10a + p10b) / (p10c + p10d) / (p10e + p10f) +
    plot_annotation(
      title = sprintf("%s -- Nuclear Temporal Trends", SPECIES),
      subtitle = sprintf("Total nuclei: %s. Depth = Edge vs Core. %d bad tps excluded. 3D%s.",
                         format(nuc_global$total_nuclei_detected, big.mark = ","),
                         length(bad_tps), core_3d_label),
      theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14),
                    plot.subtitle = element_text(hjust = 0, size = 10, color = "grey50"))
    )
  save_pdf(fig10, "zebrafish_dynamics_10_nuclear.pdf", width = 14, height = 14)

} else {
  cat("  Skipped (no nuclear stats)\n")
}

# ============================================================================
# FIGURE 4: SPATIAL CONTEXT — Zone Map & Depth Maps
# ============================================================================

cat("\n=== FIGURE 4: SPATIAL CONTEXT ===\n")

# Annotation helpers — zone boundaries
add_spatial_landmarks_zf <- function(p) {
  p <- p +
    geom_hline(yintercept = MARGIN_THETA,
               color = "gold3", linetype = "solid", linewidth = 0.8) +
    geom_hline(yintercept = ZONE_MARGIN_TOP,
               color = "gold3", linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = ZONE_PREMARG_TOP,
               color = "gold3", linetype = "dotted", linewidth = 0.4)
  p
}

# Time windows for spatial snapshots
HALF_WIN <- 15
val_centers <- round(seq(max(spots$FRAME) * 0.1, max(spots$FRAME) * 0.9, length.out = 4))

make_time_label <- function(frame, time_min) sprintf("t=%.0f min (frame %d)", time_min, frame)
time_label_levels <- make_time_label(val_centers, frame_to_min(val_centers))

facet_time <- facet_wrap(~time_label, ncol = 2)
spatial_labs <- labs(x = PHI_LABEL, y = THETA_LABEL)

# 9a: Zone map at 4 time points
val_spots <- spots %>%
  filter(FRAME %in% val_centers) %>%
  mutate(time_label = factor(make_time_label(FRAME, time_min), levels = time_label_levels))

p9a <- ggplot(val_spots, aes(x = PHI_DEG, y = THETA_DEG, color = zone)) +
  geom_point(size = 0.2, alpha = 0.4) +
  scale_color_manual(values = zone_colors, name = "Zone") +
  scale_y_reverse() + facet_time + spatial_labs +
  labs(title = "Zone assignment at four time points",
       subtitle = "Gold solid = auto-detected margin. Dashed/dotted = zone boundaries.") +
  theme_spatial(8) + theme(legend.position = "right")
p9a <- add_spatial_landmarks_zf(p9a)

# 9b: Zone composition over time
zone_time <- spots %>%
  mutate(time_bin = floor(time_min / 10) * 10) %>%
  group_by(time_bin, zone) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(time_bin) %>%
  mutate(frac = n / sum(n))

p9b <- ggplot(zone_time, aes(x = time_bin, y = frac, fill = zone)) +
  geom_area(alpha = 0.85) +
  scale_fill_manual(values = zone_colors, name = "Zone") +
  labs(title = "Zone composition over time",
       x = "Time (min)", y = "Fraction of cells") +
  theme_dynamics()

# 9c: Depth map (tile heatmap at 4 time points)
val_vel_epi <- bind_rows(lapply(val_centers, function(cf) {
  spots_vel %>%
    filter(FRAME >= (cf - HALF_WIN) & FRAME <= (cf + HALF_WIN)) %>%
    mutate(time_label = factor(make_time_label(cf, frame_to_min(cf)),
                               levels = time_label_levels))
}))

depth_tiles <- val_vel_epi %>%
  mutate(theta_bin = floor(THETA_DEG / TILE_DTHETA) * TILE_DTHETA + TILE_DTHETA / 2,
         phi_bin   = floor(PHI_DEG / TILE_DPHI) * TILE_DPHI + TILE_DPHI / 2) %>%
  group_by(time_label, theta_bin, phi_bin) %>%
  summarise(mean_depth = mean(SPHERICAL_DEPTH, na.rm = TRUE),
            v_depth = mean(v_depth, na.rm = TRUE),
            n_cells = n(), .groups = "drop")

p9c <- ggplot(depth_tiles, aes(x = phi_bin, y = theta_bin, fill = mean_depth)) +
  geom_tile(width = TILE_DPHI, height = TILE_DTHETA) +
  scale_fill_viridis(option = "inferno", name = "Depth\n(um)",
                     limits = c(0, quantile(spots$SPHERICAL_DEPTH, 0.98)), oob = scales::squish) +
  scale_y_reverse() + facet_time + spatial_labs +
  labs(title = "Spherical depth map",
       subtitle = "Blastoderm cells. Deeper = closer to yolk.") +
  theme_spatial(8) + theme(legend.position = "right")
p9c <- add_spatial_landmarks_zf(p9c)

fig9 <- (p9a) / (p9b + plot_spacer()) / (p9c) +
  plot_annotation(
    title = sprintf("%s -- Spatial Context & Zone Maps", SPECIES),
    subtitle = "Gold lines = zone boundaries.",
    theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14),
                  plot.subtitle = element_text(hjust = 0, size = 9, color = "grey40"))
  )
save_pdf(fig9, "zebrafish_dynamics_04_spatial_context.pdf", width = 16, height = 20)

# ============================================================================
# FIGURE 5: SPATIAL FLOW SUMMARY — Vector Fields & Tile Maps
# ============================================================================

cat("\n=== FIGURE 5: SPATIAL FLOW SUMMARY ===\n")

flow_field <- spots_vel %>%
  mutate(theta_bin = floor(THETA_DEG / THETA_BIN) * THETA_BIN + THETA_BIN / 2,
         phi_bin   = floor(PHI_DEG / PHI_BIN) * PHI_BIN + PHI_BIN / 2) %>%
  group_by(theta_bin, phi_bin) %>%
  summarise(
    mean_v_epi   = mean(v_epiboly, na.rm = TRUE),
    mean_v_lat   = mean(v_convergence, na.rm = TRUE),
    mean_v_depth = mean(v_depth, na.rm = TRUE),
    mean_speed   = mean(inst_speed, na.rm = TRUE),
    n = n(), .groups = "drop"
  ) %>%
  filter(n >= FLOW_MIN_N)

# 10a: Velocity vector field
p10a <- ggplot(flow_field) +
  geom_segment(aes(x = phi_bin, y = theta_bin,
                   xend = phi_bin + mean_v_lat * 3,
                   yend = theta_bin + mean_v_epi * 3,
                   color = mean_speed),
               arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
               linewidth = 0.4) +
  scale_color_viridis(option = "B", name = "Speed\n(um/min)") +
  scale_y_reverse() + spatial_labs +
  labs(title = "Velocity vector field",
       subtitle = "Arrow = tangential velocity (epiboly + lateral). Color = speed.") +
  theme_spatial()
p10a <- add_spatial_landmarks_zf(p10a)

# 10b: Speed map
p10b <- ggplot(flow_field, aes(x = phi_bin, y = theta_bin, fill = mean_speed)) +
  geom_tile(width = PHI_BIN, height = THETA_BIN) +
  scale_fill_viridis(option = "B", name = "Speed\n(um/min)") +
  scale_y_reverse() + spatial_labs +
  labs(title = "Mean speed on embryo surface") +
  theme_spatial()
p10b <- add_spatial_landmarks_zf(p10b)

# 10c: Depth velocity map
p10c <- ggplot(flow_field, aes(x = phi_bin, y = theta_bin, fill = mean_v_depth)) +
  geom_tile(width = PHI_BIN, height = THETA_BIN) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey20", high = "#D73027",
                       midpoint = 0, limits = c(-0.05, 0.05), oob = scales::squish,
                       name = "V_depth\n(um/min)") +
  scale_y_reverse() + spatial_labs +
  labs(title = "Depth change rate",
       subtitle = "Red = deepening. Blue = surfacing.") +
  theme_spatial()
p10c <- add_spatial_landmarks_zf(p10c)

# 10d: Epiboly velocity map
p10d <- ggplot(flow_field, aes(x = phi_bin, y = theta_bin, fill = mean_v_epi)) +
  geom_tile(width = PHI_BIN, height = THETA_BIN) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey20", high = "#EF8A62",
                       midpoint = 0, name = "V_epiboly\n(um/min)") +
  scale_y_reverse() + spatial_labs +
  labs(title = "Epiboly velocity",
       subtitle = "Orange = vegetal-ward. Blue = animal-ward.") +
  theme_spatial()
p10d <- add_spatial_landmarks_zf(p10d)

# 10e: Movement type spatial map
type_spatial <- spots_vel %>%
  left_join(track_metrics %>% dplyr::select(TRACK_ID, movement_type), by = "TRACK_ID") %>%
  filter(!is.na(movement_type)) %>%
  mutate(theta_bin = floor(THETA_DEG / THETA_BIN) * THETA_BIN + THETA_BIN / 2,
         phi_bin   = floor(PHI_DEG / PHI_BIN) * PHI_BIN + PHI_BIN / 2) %>%
  group_by(theta_bin, phi_bin) %>%
  summarise(
    n_epi  = sum(movement_type == "Epiboly"),
    n_ani  = sum(movement_type == "Animalward"),
    n_conv = sum(movement_type == "Convergence"),
    n_ingr = sum(movement_type == "Ingression"),
    n_total = n(), .groups = "drop"
  ) %>%
  filter(n_total >= FLOW_MIN_N) %>%
  mutate(dominant = case_when(
    n_ingr >= n_epi & n_ingr >= n_ani & n_ingr >= n_conv ~ "Ingression",
    n_conv >= n_epi & n_conv >= n_ani                    ~ "Convergence",
    n_epi >= n_ani                                       ~ "Epiboly",
    TRUE                                                 ~ "Animalward"
  ))

# Also compute per-tile mean velocity for arrow overlay
type_arrows <- spots_vel %>%
  mutate(theta_bin = floor(THETA_DEG / THETA_BIN) * THETA_BIN + THETA_BIN / 2,
         phi_bin   = floor(PHI_DEG / PHI_BIN) * PHI_BIN + PHI_BIN / 2) %>%
  group_by(theta_bin, phi_bin) %>%
  summarise(mean_v_epi = mean(v_epiboly, na.rm = TRUE),
            mean_v_lat = mean(v_convergence, na.rm = TRUE),
            mean_speed = mean(inst_speed, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= FLOW_MIN_N)

p10e <- ggplot(type_spatial, aes(x = phi_bin, y = theta_bin, fill = dominant)) +
  geom_tile(width = PHI_BIN, height = THETA_BIN) +
  geom_segment(data = type_arrows,
               aes(x = phi_bin, y = theta_bin,
                   xend = phi_bin + mean_v_lat * 3,
                   yend = theta_bin + mean_v_epi * 3),
               arrow = arrow(length = unit(0.06, "cm"), type = "closed"),
               linewidth = 0.3, color = "white", alpha = 0.7,
               inherit.aes = FALSE) +
  scale_fill_manual(values = mvt_colors, name = "Dominant\ntype") +
  scale_y_reverse() + spatial_labs +
  labs(title = "Dominant movement type per tile",
       subtitle = "Arrows = mean velocity direction.") +
  theme_spatial()
p10e <- add_spatial_landmarks_zf(p10e)

# 10f: Flow field at 3 time points
time_points_ff <- c(round(max(spots$FRAME) * 0.1),
                    round(max(spots$FRAME) * 0.5),
                    round(max(spots$FRAME) * 0.9))
time_labels_ordered <- sprintf("t = %.0f min", frame_to_min(time_points_ff))

flow_timepoints <- lapply(seq_along(time_points_ff), function(i) {
  t <- time_points_ff[i]
  spots_vel %>%
    filter(FRAME >= t - 5 & FRAME <= t + 5) %>%
    mutate(theta_bin = floor(THETA_DEG / (THETA_BIN * 2)) * (THETA_BIN * 2) + THETA_BIN,
           phi_bin   = floor(PHI_DEG / (PHI_BIN * 2)) * (PHI_BIN * 2) + PHI_BIN) %>%
    group_by(theta_bin, phi_bin) %>%
    summarise(mean_v_epi = mean(v_epiboly, na.rm = TRUE),
              mean_v_lat = mean(v_convergence, na.rm = TRUE),
              mean_speed = mean(inst_speed, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    filter(n >= 5) %>%
    mutate(time_label = factor(time_labels_ordered[i], levels = time_labels_ordered))
}) %>% bind_rows()

p10f <- ggplot(flow_timepoints) +
  geom_segment(aes(x = phi_bin, y = theta_bin,
                   xend = phi_bin + mean_v_lat * 3,
                   yend = theta_bin + mean_v_epi * 3,
                   color = mean_speed),
               arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
               linewidth = 0.4) +
  scale_color_viridis(option = "B", name = "Speed") +
  facet_wrap(~time_label) + scale_y_reverse() + spatial_labs +
  labs(title = "Flow fields at three stages") +
  theme_spatial()

# 10g: Movement type evolution over time
type_time <- spots_vel %>%
  left_join(track_metrics %>% dplyr::select(TRACK_ID, movement_type), by = "TRACK_ID") %>%
  filter(!is.na(movement_type)) %>%
  mutate(time_bin = floor(time_min / 5) * 5) %>%
  group_by(time_bin, movement_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(time_bin) %>%
  mutate(frac = n / sum(n))

p10g <- ggplot(type_time, aes(x = time_bin, y = frac, fill = movement_type)) +
  geom_area() +
  scale_fill_manual(values = mvt_colors, name = "Type") +
  labs(title = "Movement type fraction over time", x = "Time (min)", y = "Fraction") +
  theme_dynamics()

fig10 <- (p10a + p10b) / (p10c + p10d) / (p10e + p10g) / (p10f) +
  plot_annotation(
    title = sprintf("%s -- Spatial Flow Summary", SPECIES),
    subtitle = "Time-averaged velocity fields on embryo surface. Gold = zone boundaries.",
    theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14),
                  plot.subtitle = element_text(hjust = 0, size = 8, color = "grey40"))
  )
save_pdf(fig10, "zebrafish_dynamics_05_spatial_flow.pdf", width = 14, height = 20)

# =============================================================================
# SAVE SUMMARY
# =============================================================================

cat("\n=== SAVING SUMMARY ===\n")

write_csv(track_metrics, file.path(OUTPUT_DIR, "zebrafish_dynamics_track_metrics.csv"))

summary_data <- data.frame(
  species = SPECIES,
  metric = c("n_spots", "n_tracks", "n_frames", "frame_interval_sec",
             "mean_speed", "median_speed", "q95_speed",
             "mean_speed_moving", "median_speed_moving", "q95_speed_moving",
             "pct_zero_disp",
             "mean_v_epiboly", "mean_v_convergence", "mean_v_depth",
             "msd_alpha", "mean_straightness", "median_straightness",
             "mean_turning_angle", "sphere_radius_um",
             "mean_nuclear_diam_um"),
  value = c(nrow(spots), length(unique(spots$TRACK_ID)),
            max(spots$FRAME) - min(spots$FRAME) + 1, FRAME_INTERVAL_SEC,
            speed_summary$mean_speed, speed_summary$median_speed, speed_summary$q95_speed,
            speed_summary_moving$mean_speed, speed_summary_moving$median_speed, speed_summary_moving$q95_speed,
            100 * n_zero_disp / nrow(spots_vel),
            mean(spots_vel$v_epiboly, na.rm = TRUE),
            mean(spots_vel$v_convergence, na.rm = TRUE),
            mean(spots_vel$v_depth, na.rm = TRUE),
            alpha_global,
            mean(track_metrics$straightness, na.rm = TRUE),
            median(track_metrics$straightness, na.rm = TRUE),
            mean(track_metrics$mean_turning, na.rm = TRUE),
            SPHERE_RADIUS,
            ifelse(HAS_NUCLEAR, mean(nuc_time$eqsph_diam_mean_um, na.rm = TRUE), NA_real_))
)
write_csv(summary_data, file.path(OUTPUT_DIR, "zebrafish_dynamics_summary.csv"))

# Save thickness for comparison
write_csv(thickness_zone_time %>% mutate(species = SPECIES),
          file.path(OUTPUT_DIR, "zebrafish_thickness_data.csv"))

cat("\n")
cat(strrep("=", 70), "\n")
cat("  ZEBRAFISH ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n")
