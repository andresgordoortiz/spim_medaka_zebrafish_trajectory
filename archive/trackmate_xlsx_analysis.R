# =============================================================================
# TrackMate Analysis from Imaris CSV Export
# =============================================================================
#
# This script reads tracking data from Imaris CSV exports in the all_tracks/
# folder and performs the same analyses as trackmate_analysis.R.
#
# DATA STRUCTURE:
#   Imaris exports each statistic as a separate CSV file.
#   Each CSV has 4 header lines:
#     Line 1: Empty
#     Line 2: Statistic name (e.g., "Position")
#     Line 3: Separator (" ==================== ")
#     Line 4: Column names
#     Lines 5+: Data
#
#   Files (all prefixed ome-tiff.companion-bin-8bit-crop_):
#     Position.csv                → X, Y, Z, Time (frame), TrackID, SpotID
#     Track_Duration.csv          → duration per track (IN FRAMES, not seconds!)
#     Track_Length.csv             → total path length per track (µm)
#     Track_Speed_Mean.csv        → mean speed (µm/frame, NOT µm/s!)
#     Track_Straightness.csv      → confinement ratio (unitless)
#     Track_Displacement_Length.csv → net displacement per track (µm)
#     Track_Displacement.csv      → displacement X/Y/Z (redundant)
#     Track_Displacement_X/Y/Z.csv → individual components (redundant)
#
#   CRITICAL: Imaris uses UNCALIBRATED time (1 frame = 1 second).
#   We correct all time-dependent values by FRAME_INTERVAL_SEC.
#
# =============================================================================

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(viridis)
library(purrr)
library(tibble)
library(patchwork)
library(cluster)
library(plotly)
library(htmlwidgets)

# =============================================================================
# GLOBAL PARAMETERS
# =============================================================================

DATA_DIR <- "all_tracks"
FILE_PREFIX <- "ome-tiff.companion-bin-8bit-crop"
FRAME_INTERVAL_SEC <- 120  # Actual seconds per frame — CHANGE THIS to match your acquisition
FRAME_INTERVAL_MIN <- FRAME_INTERVAL_SEC / 60
OUTPUT_DIR <- "analysis_output_imaris"

# =============================================================================
# PLOT STYLING
# =============================================================================

movement_colors <- c(
  "Directed" = "#2166AC",
  "Confined" = "#B2182B",
  "Random"   = "#762A83"
)

theme_track <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 1),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1, color = "grey50"),
      legend.position = "bottom",
      legend.key.size = unit(0.3, "cm"),
      legend.text = element_text(size = base_size - 2),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = base_size - 1)
    )
}

# =============================================================================
# PART 1: DATA LOADING
# =============================================================================
#
# Each Imaris CSV has 4 header lines:
#   Line 1: empty
#   Line 2: Statistic name (e.g., "Position")
#   Line 3: Separator (" ==================== ")
#   Line 4: Column names
# We skip 3 lines so line 4 becomes the header.
# Trailing commas in each row produce an extra empty column — we drop it.
#

cat("\n=== PART 1: DATA LOADING FROM IMARIS CSVs ===\n")
cat(sprintf("  Folder: %s/\n", DATA_DIR))
cat(sprintf("  Frame interval: %d sec (%.1f min)\n", FRAME_INTERVAL_SEC, FRAME_INTERVAL_MIN))

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# Helper: read an Imaris CSV (skip 3 header lines, drop trailing empty column)
read_imaris_csv <- function(stat_name) {
  path <- file.path(DATA_DIR, paste0(FILE_PREFIX, "_", stat_name, ".csv"))
  df <- read_csv(path, skip = 3, show_col_types = FALSE)
  # Drop trailing empty column from trailing comma
  df <- df %>% select(where(~!all(is.na(.x))))
  df
}

# --- Spots data ---
spots <- read_imaris_csv("Position") %>%
  rename(
    POSITION_X = `Position X`,
    POSITION_Y = `Position Y`,
    POSITION_Z = `Position Z`,
    FRAME      = Time,
    TRACK_ID   = TrackID,
    ID         = ID
  ) %>%
  select(POSITION_X, POSITION_Y, POSITION_Z, FRAME, TRACK_ID, ID) %>%
  mutate(across(everything(), as.numeric))

# --- Track-level data ---
track_duration <- read_imaris_csv("Track_Duration") %>%
  rename(TRACK_DURATION_FRAMES = `Track Duration`, TRACK_ID = ID) %>%
  select(TRACK_ID, TRACK_DURATION_FRAMES) %>%
  mutate(across(everything(), as.numeric))

track_length <- read_imaris_csv("Track_Length") %>%
  rename(TRACK_LENGTH = `Track Length`, TRACK_ID = ID) %>%
  select(TRACK_ID, TRACK_LENGTH) %>%
  mutate(across(everything(), as.numeric))

track_speed <- read_imaris_csv("Track_Speed_Mean") %>%
  rename(TRACK_SPEED_MEAN_RAW = `Track Speed Mean`, TRACK_ID = ID) %>%
  select(TRACK_ID, TRACK_SPEED_MEAN_RAW) %>%
  mutate(across(everything(), as.numeric))

track_straightness <- read_imaris_csv("Track_Straightness") %>%
  rename(CONFINEMENT_RATIO = `Track Straightness`, TRACK_ID = ID) %>%
  select(TRACK_ID, CONFINEMENT_RATIO) %>%
  mutate(across(everything(), as.numeric))

track_displacement <- read_imaris_csv("Track_Displacement_Length") %>%
  rename(TRACK_DISPLACEMENT = `Track Displacement Length`, TRACK_ID = ID) %>%
  select(TRACK_ID, TRACK_DISPLACEMENT) %>%
  mutate(across(everything(), as.numeric))

# Merge all track-level data
tracks <- track_duration %>%
  left_join(track_length,       by = "TRACK_ID") %>%
  left_join(track_speed,        by = "TRACK_ID") %>%
  left_join(track_straightness, by = "TRACK_ID") %>%
  left_join(track_displacement, by = "TRACK_ID")

# Convert units
# IMPORTANT: Imaris uses UNCALIBRATED time (1 frame = 1 "second").
#   "Track Duration" value is actually N_FRAMES (not real seconds)
#   "Track Speed Mean" is path_length / N_FRAMES (µm/frame, not µm/s)
# Correction:
#   Real duration (min) = N_FRAMES * FRAME_INTERVAL_SEC / 60
#   Real speed (µm/min) = (µm/frame) / FRAME_INTERVAL_SEC * 60
tracks <- tracks %>%
  mutate(
    DURATION_MIN      = TRACK_DURATION_FRAMES * FRAME_INTERVAL_SEC / 60,
    SPEED_MEAN_UM_MIN = TRACK_SPEED_MEAN_RAW / FRAME_INTERVAL_SEC * 60
  )

# Count spots per track
spots_per_track <- spots %>% count(TRACK_ID, name = "N_SPOTS")
tracks <- tracks %>% left_join(spots_per_track, by = "TRACK_ID")

cat(sprintf("  Spots:  %d detections\n", nrow(spots)))
cat(sprintf("  Tracks: %d\n", nrow(tracks)))
cat(sprintf("  Frames: %d (%.1f min total)\n",
            max(spots$FRAME), max(spots$FRAME) * FRAME_INTERVAL_MIN))
cat(sprintf("  Frame interval: %d sec\n\n", FRAME_INTERVAL_SEC))

# =============================================================================
# PART 2: DATA EXPLORATION
# =============================================================================

cat("=== PART 2: DATA EXPLORATION ===\n\n")

track_summary <- tracks %>%
  summarise(
    duration_min_median  = median(DURATION_MIN, na.rm = TRUE),
    speed_um_min_median  = median(SPEED_MEAN_UM_MIN, na.rm = TRUE),
    displacement_median  = median(TRACK_DISPLACEMENT, na.rm = TRUE),
    confinement_median   = median(CONFINEMENT_RATIO, na.rm = TRUE),
    n_spots_median       = median(N_SPOTS, na.rm = TRUE),
    track_length_median  = median(TRACK_LENGTH, na.rm = TRUE)
  )

cat("Track Statistics:\n")
cat(sprintf("  Median duration: %.1f min\n", track_summary$duration_min_median))
cat(sprintf("  Median speed: %.2f µm/min\n", track_summary$speed_um_min_median))
cat(sprintf("  Median displacement: %.1f µm\n", track_summary$displacement_median))
cat(sprintf("  Median confinement ratio: %.2f\n", track_summary$confinement_median))
cat(sprintf("  Median track length: %.1f µm\n", track_summary$track_length_median))
cat(sprintf("  Median n spots: %.0f\n\n", track_summary$n_spots_median))

# =============================================================================
# PART 3: QUALITY CONTROL
# =============================================================================
#
# Filter criteria:
#   - Remove very short tracks (<3 spots)
#   - Remove unrealistically fast tracks (mean speed >10 µm/min)
#

cat("=== PART 3: QUALITY CONTROL ===\n\n")

MIN_SPOTS <- 3
MAX_SPEED_UM_MIN <- 10

tracks <- tracks %>%
  mutate(
    flag_too_short = N_SPOTS < MIN_SPOTS,
    flag_too_fast  = SPEED_MEAN_UM_MIN > MAX_SPEED_UM_MIN,
    passes_qc      = !flag_too_short & !flag_too_fast
  )

cat(sprintf("  Too short (<%d spots): %d tracks\n", MIN_SPOTS, sum(tracks$flag_too_short)))
cat(sprintf("  Too fast (>%d µm/min mean): %d tracks\n", MAX_SPEED_UM_MIN, sum(tracks$flag_too_fast)))

tracks_clean <- tracks %>% filter(passes_qc)
spots_clean  <- spots %>% filter(TRACK_ID %in% tracks_clean$TRACK_ID)

cat(sprintf("  Kept: %d / %d tracks (%.1f%%)\n",
            nrow(tracks_clean), nrow(tracks), 100 * nrow(tracks_clean) / nrow(tracks)))
cat(sprintf("  Kept: %d / %d spots\n\n", nrow(spots_clean), nrow(spots)))

# =============================================================================
# PART 4: QC VISUALIZATIONS
# =============================================================================

cat("=== PART 4: QC VISUALIZATIONS ===\n")

# Track length distribution
p1 <- ggplot(tracks, aes(x = N_SPOTS, fill = passes_qc)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "grey70"),
                    labels = c("TRUE" = "Kept", "FALSE" = "Filtered")) +
  labs(title = "Track Length Distribution",
       subtitle = sprintf("Min %d spots filter | %d / %d tracks kept", MIN_SPOTS,
                          nrow(tracks_clean), nrow(tracks)),
       x = "Number of spots", y = "Count", fill = NULL) +
  theme_track()

# Speed distribution
p2 <- ggplot(tracks_clean, aes(x = SPEED_MEAN_UM_MIN)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "#2166AC", alpha = 0.7) +
  geom_density(color = "#B2182B", linewidth = 1) +
  labs(title = "Mean Track Speed",
       subtitle = sprintf("Median: %.2f µm/min", median(tracks_clean$SPEED_MEAN_UM_MIN)),
       x = "Mean speed (µm/min)", y = "Density") +
  theme_track()

# Duration vs speed
p3 <- ggplot(tracks, aes(x = DURATION_MIN, y = SPEED_MEAN_UM_MIN, color = passes_qc)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                     labels = c("TRUE" = "Kept", "FALSE" = "Filtered")) +
  labs(title = "Duration vs Speed",
       subtitle = "Identify outliers",
       x = "Duration (min)", y = "Mean speed (µm/min)", color = NULL) +
  theme_track()

# Confinement distribution raw
p4 <- ggplot(tracks_clean, aes(x = CONFINEMENT_RATIO)) +
  geom_histogram(bins = 25, fill = "#2166AC", alpha = 0.7) +
  labs(title = "Track Straightness (raw)",
       x = "Straightness (confinement ratio)", y = "Count") +
  theme_track()

qc_combined <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Quality Control",
    subtitle = sprintf("Imaris data: %d tracks, %d spots", nrow(tracks), nrow(spots)),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey50")
    )
  )

ggsave(file.path(OUTPUT_DIR, "01_quality_control.pdf"), qc_combined, width = 10, height = 10)
cat("  Saved: 01_quality_control.pdf\n")

# =============================================================================
# PART 5: VELOCITY ANALYSIS
# =============================================================================

cat("\n=== PART 5: VELOCITY ANALYSIS ===\n")

spots_velocity <- spots_clean %>%
  arrange(TRACK_ID, FRAME) %>%
  group_by(TRACK_ID) %>%
  mutate(
    dt_min       = (FRAME - lag(FRAME)) * FRAME_INTERVAL_MIN,
    dx           = POSITION_X - lag(POSITION_X),
    dy           = POSITION_Y - lag(POSITION_Y),
    dz           = POSITION_Z - lag(POSITION_Z),
    displacement = sqrt(dx^2 + dy^2 + dz^2),
    inst_speed_um_min = displacement / dt_min,
    vx = dx / dt_min,
    vy = dy / dt_min,
    vz = dz / dt_min,
    time_min = (FRAME - min(FRAME)) * FRAME_INTERVAL_MIN
  ) %>%
  ungroup()

speed_stats <- spots_velocity %>%
  filter(!is.na(inst_speed_um_min)) %>%
  summarise(
    mean   = mean(inst_speed_um_min),
    median = median(inst_speed_um_min),
    sd     = sd(inst_speed_um_min),
    q95    = quantile(inst_speed_um_min, 0.95)
  )

cat(sprintf("  Mean: %.2f µm/min | Median: %.2f | 95th pct: %.2f\n",
            speed_stats$mean, speed_stats$median, speed_stats$q95))

# Instantaneous speed distribution
p5 <- ggplot(spots_velocity %>% filter(!is.na(inst_speed_um_min)),
             aes(x = inst_speed_um_min)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "#2166AC", alpha = 0.7) +
  geom_density(color = "#B2182B", linewidth = 1) +
  geom_vline(xintercept = speed_stats$median, linetype = "dashed", color = "grey30", linewidth = 0.4) +
  labs(
    title = "Instantaneous Nuclear Speed (3D)",
    subtitle = sprintf("Median: %.2f µm/min | 95th pct: %.2f µm/min",
                       speed_stats$median, speed_stats$q95),
    x = "Speed (µm/min)", y = "Density"
  ) +
  theme_track()

# Speed over time
speed_over_time <- spots_velocity %>%
  filter(!is.na(inst_speed_um_min)) %>%
  mutate(time_bin_min = floor(FRAME * FRAME_INTERVAL_MIN / 10) * 10) %>%
  group_by(time_bin_min) %>%
  summarise(
    mean_speed = mean(inst_speed_um_min, na.rm = TRUE),
    se_speed   = sd(inst_speed_um_min, na.rm = TRUE) / sqrt(n()),
    n = n()
  )

p6 <- ggplot(speed_over_time, aes(x = time_bin_min, y = mean_speed)) +
  geom_ribbon(aes(ymin = mean_speed - se_speed, ymax = mean_speed + se_speed),
              alpha = 0.3, fill = "#2166AC") +
  geom_line(color = "#2166AC", linewidth = 1) +
  labs(
    title = "Mean Instantaneous Speed Over Time",
    subtitle = "Shaded: SE",
    x = "Time (min from start)", y = "Mean speed (µm/min)"
  ) +
  theme_track()

velocity_combined <- (p5 + p6) +
  plot_annotation(
    title = "Velocity Analysis",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "02_velocity_analysis.pdf"), velocity_combined, width = 14, height = 5)
cat("  Saved: 02_velocity_analysis.pdf\n")

# =============================================================================
# PART 6: DIRECTIONALITY ANALYSIS
# =============================================================================

cat("\n=== PART 6: DIRECTIONALITY ANALYSIS ===\n")

spots_direction <- spots_velocity %>%
  filter(!is.na(dx) & !is.na(dy) & !is.na(dz)) %>%
  group_by(TRACK_ID) %>%
  mutate(
    prev_dx = lag(dx), prev_dy = lag(dy), prev_dz = lag(dz),
    dot_product  = dx * prev_dx + dy * prev_dy + dz * prev_dz,
    mag_current  = sqrt(dx^2 + dy^2 + dz^2),
    mag_prev     = sqrt(prev_dx^2 + prev_dy^2 + prev_dz^2),
    turning_angle = acos(pmin(1, pmax(-1, dot_product / (mag_current * mag_prev + 1e-10)))) * 180 / pi,
    direction_xy  = atan2(dy, dx) * 180 / pi
  ) %>%
  ungroup()

turning_stats <- spots_direction %>%
  filter(!is.na(turning_angle)) %>%
  summarise(mean_angle = mean(turning_angle), median_angle = median(turning_angle))

cat(sprintf("  Mean turning angle: %.1f° | Median: %.1f°\n",
            turning_stats$mean_angle, turning_stats$median_angle))

p7 <- ggplot(spots_direction %>% filter(!is.na(turning_angle)),
             aes(x = turning_angle)) +
  geom_histogram(aes(y = after_stat(density)), bins = 36, fill = "#2166AC", alpha = 0.7) +
  geom_vline(xintercept = 90, linetype = "dashed", color = "#B2182B", linewidth = 0.4) +
  labs(title = "Turning Angles", subtitle = "90° = random; <90° = persistent; >90° = reversing",
       x = "Turning angle (°)", y = "Density") +
  theme_track()

p8 <- ggplot(spots_direction %>% filter(!is.na(direction_xy)),
             aes(x = direction_xy)) +
  geom_histogram(bins = 36, fill = "#2166AC", color = "white", alpha = 0.7) +
  coord_polar(start = 0) +
  scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 135, 45)) +
  labs(title = "Movement Direction (XY)", subtitle = "Rose diagram",
       x = "Direction (°)", y = "Count") +
  theme_track() +
  theme(axis.text.y = element_blank())

direction_combined <- (p7 + p8) +
  plot_annotation(
    title = "Directionality Analysis",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "03_directionality.pdf"), direction_combined, width = 12, height = 5)
cat("  Saved: 03_directionality.pdf\n")

# =============================================================================
# PART 7: SPATIAL VELOCITY PATTERNS
# =============================================================================

cat("\n=== PART 7: SPATIAL VELOCITY PATTERNS ===\n")

BIN_SIZE <- 30  # µm

spatial_velocity <- spots_velocity %>%
  filter(!is.na(inst_speed_um_min)) %>%
  mutate(
    x_bin = floor(POSITION_X / BIN_SIZE) * BIN_SIZE,
    y_bin = floor(POSITION_Y / BIN_SIZE) * BIN_SIZE
  )

local_stats_xy <- spatial_velocity %>%
  group_by(x_bin, y_bin) %>%
  summarise(
    mean_speed = mean(inst_speed_um_min, na.rm = TRUE),
    mean_vx    = mean(vx, na.rm = TRUE),
    mean_vy    = mean(vy, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

p9 <- ggplot(local_stats_xy, aes(x = x_bin, y = y_bin, fill = mean_speed)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Speed\n(µm/min)", option = "inferno") +
  coord_fixed() +
  labs(title = "Spatial Speed Heatmap (XY)", x = "X (µm)", y = "Y (µm)") +
  theme_track()

p10 <- ggplot(local_stats_xy %>% filter(n >= 10),
              aes(x = x_bin, y = y_bin)) +
  geom_segment(aes(xend = x_bin + mean_vx * 5,
                   yend = y_bin + mean_vy * 5,
                   color = mean_speed),
               arrow = arrow(length = unit(0.08, "cm")), linewidth = 0.4) +
  scale_color_viridis_c(name = "Speed", option = "inferno") +
  coord_fixed() +
  labs(title = "Velocity Vector Field (XY)", x = "X (µm)", y = "Y (µm)") +
  theme_track()

spatial_combined <- (p9 + p10) +
  plot_annotation(
    title = "Spatial Velocity Patterns",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "04_spatial_velocity.pdf"), spatial_combined, width = 12, height = 6)
cat("  Saved: 04_spatial_velocity.pdf\n")

# =============================================================================
# PART 8: MSD ANALYSIS
# =============================================================================

cat("\n=== PART 8: MSD ANALYSIS ===\n")

calculate_msd <- function(df, max_lag = 30) {
  df <- df %>% arrange(FRAME)
  n <- nrow(df)
  if (n < 3) return(NULL)
  max_lag <- min(max_lag, n - 1)
  msd_values <- sapply(1:max_lag, function(lag) {
    displacements <- (df$POSITION_X[(lag+1):n] - df$POSITION_X[1:(n-lag)])^2 +
      (df$POSITION_Y[(lag+1):n] - df$POSITION_Y[1:(n-lag)])^2 +
      (df$POSITION_Z[(lag+1):n] - df$POSITION_Z[1:(n-lag)])^2
    mean(displacements, na.rm = TRUE)
  })
  tibble(lag_frames = 1:max_lag, lag_min = 1:max_lag * FRAME_INTERVAL_MIN, msd = msd_values)
}

# Sample tracks for MSD (computing MSD for all 5000+ tracks is very slow)
MAX_TRACKS_MSD <- 500
msd_track_ids <- unique(spots_clean$TRACK_ID)
if (length(msd_track_ids) > MAX_TRACKS_MSD) {
  set.seed(42)
  msd_track_ids <- sample(msd_track_ids, MAX_TRACKS_MSD)
  cat(sprintf("  Sampling %d / %d tracks for MSD analysis\n",
              MAX_TRACKS_MSD, n_distinct(spots_clean$TRACK_ID)))
}

msd_data <- spots_clean %>%
  filter(TRACK_ID %in% msd_track_ids) %>%
  arrange(TRACK_ID, FRAME) %>%
  group_by(TRACK_ID) %>%
  group_modify(~{
    result <- calculate_msd(.x)
    if (is.null(result)) return(tibble())
    result
  }) %>%
  ungroup()

ensemble_msd <- msd_data %>%
  group_by(lag_min) %>%
  summarise(
    mean_msd = mean(msd, na.rm = TRUE),
    se_msd   = sd(msd, na.rm = TRUE) / sqrt(n()),
    n = n()
  ) %>%
  filter(n >= 10)

# Fit power law
msd_fit <- tryCatch({
  lm(log(mean_msd) ~ log(lag_min), data = ensemble_msd %>% filter(lag_min <= 10))
}, error = function(e) NULL)

if (!is.null(msd_fit)) {
  alpha <- coef(msd_fit)[2]
  D_apparent <- exp(coef(msd_fit)[1]) / 6
  cat(sprintf("  α = %.3f → %s\n", alpha,
              ifelse(alpha < 0.8, "subdiffusive", ifelse(alpha > 1.2, "superdiffusive", "diffusive"))))
  cat(sprintf("  D_apparent = %.3f µm²/min\n", D_apparent))
} else {
  alpha <- NA
  D_apparent <- NA
  cat("  Could not fit MSD power law\n")
}

# MSD plot
p11 <- ggplot() +
  geom_line(data = msd_data %>%
              filter(TRACK_ID %in% sample(unique(TRACK_ID), min(50, n_distinct(TRACK_ID)))),
            aes(x = lag_min, y = msd, group = TRACK_ID), alpha = 0.2, color = "grey50") +
  geom_ribbon(data = ensemble_msd,
              aes(x = lag_min, ymin = mean_msd - se_msd, ymax = mean_msd + se_msd),
              alpha = 0.3, fill = "#B2182B") +
  geom_line(data = ensemble_msd, aes(x = lag_min, y = mean_msd),
            color = "#B2182B", linewidth = 1.2) +
  scale_x_log10() + scale_y_log10() +
  labs(
    title = "Mean Squared Displacement Analysis",
    subtitle = if (!is.na(alpha)) sprintf("α = %.2f (%s)",
                  alpha, ifelse(alpha < 0.8, "subdiffusive",
                                ifelse(alpha > 1.2, "superdiffusive", "diffusive"))) else "Fit failed",
    x = "Time lag (min, log)", y = "MSD (µm², log)"
  ) +
  theme_track()

ggsave(file.path(OUTPUT_DIR, "05_msd_analysis.pdf"), p11, width = 8, height = 6)
cat("  Saved: 05_msd_analysis.pdf\n")

# =============================================================================
# PART 9: MOVEMENT TYPE CLASSIFICATION
# =============================================================================

cat("\n=== PART 9: MOVEMENT TYPE CLASSIFICATION ===\n")

# Note: Imaris provides "Track Straightness" which = confinement ratio
# We don't have LINEARITY_OF_FORWARD_PROGRESSION from Imaris,
# so we classify based on confinement ratio + speed instead

tracks_classified <- tracks_clean %>%
  mutate(
    movement_type = case_when(
      CONFINEMENT_RATIO > 0.5 & SPEED_MEAN_UM_MIN > median(SPEED_MEAN_UM_MIN) ~ "Directed",
      CONFINEMENT_RATIO < 0.25 ~ "Confined",
      TRUE ~ "Random"
    ),
    movement_type = factor(movement_type, levels = c("Directed", "Random", "Confined"))
  )

movement_summary <- tracks_classified %>%
  count(movement_type) %>%
  mutate(pct = 100 * n / sum(n))

cat("Movement Type Classification:\n")
for (i in 1:nrow(movement_summary)) {
  cat(sprintf("  %s: %d tracks (%.1f%%)\n",
              movement_summary$movement_type[i], movement_summary$n[i], movement_summary$pct[i]))
}

# Scatter: confinement vs speed (replacing confinement vs linearity)
p12 <- ggplot(tracks_classified,
              aes(x = CONFINEMENT_RATIO, y = SPEED_MEAN_UM_MIN, color = movement_type)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_vline(xintercept = c(0.25, 0.5), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = movement_colors) +
  labs(title = "Movement Classification",
       subtitle = "Based on straightness & speed",
       x = "Straightness (confinement ratio)", y = "Mean speed (µm/min)", color = NULL) +
  theme_track()

p13 <- ggplot(movement_summary, aes(x = movement_type, y = pct, fill = movement_type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)), vjust = -0.2, size = 3, fontface = "bold") +
  scale_fill_manual(values = movement_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Movement Type Distribution", x = NULL, y = "% of tracks", fill = NULL) +
  theme_track() +
  theme(legend.position = "none", panel.grid.major.x = element_blank())

classification_combined <- (p12 + p13) +
  plot_annotation(
    title = "Nuclear Movement Classification",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "06_movement_classification.pdf"), classification_combined, width = 12, height = 5)
cat("  Saved: 06_movement_classification.pdf\n")

# =============================================================================
# PART 10: CONFINEMENT ANALYSIS
# =============================================================================

cat("\n=== PART 10: CONFINEMENT ANALYSIS ===\n")

p15 <- ggplot(tracks_clean, aes(x = CONFINEMENT_RATIO)) +
  geom_histogram(bins = 25, fill = "#2166AC", alpha = 0.7) +
  geom_vline(xintercept = c(0.3, 0.7), linetype = "dashed", color = "#B2182B", linewidth = 0.4) +
  labs(title = "Confinement Ratio Distribution",
       subtitle = "<0.3 = confined; >0.7 = directed; between = random",
       x = "Confinement ratio (straightness)", y = "Count") +
  theme_track()

p16 <- ggplot(tracks_clean, aes(x = CONFINEMENT_RATIO, y = SPEED_MEAN_UM_MIN)) +
  geom_point(alpha = 0.5, size = 2, color = "#2166AC") +
  geom_smooth(method = "loess", color = "#B2182B", linewidth = 1, se = FALSE) +
  labs(title = "Confinement vs Speed",
       subtitle = "Are faster nuclei more directed?",
       x = "Confinement ratio", y = "Mean speed (µm/min)") +
  theme_track()

confinement_combined <- (p15 + p16) +
  plot_annotation(
    title = "Confinement Analysis",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "08_confinement.pdf"), confinement_combined, width = 12, height = 5)
cat("  Saved: 08_confinement.pdf\n")

# =============================================================================
# PART 11: Z-DEPTH ANALYSIS
# =============================================================================

cat("\n=== PART 11: Z-DEPTH ANALYSIS ===\n")

z_analysis <- spots_velocity %>%
  filter(!is.na(inst_speed_um_min)) %>%
  mutate(z_layer = cut(POSITION_Z, breaks = 5,
                       labels = c("Deep", "Mid-deep", "Middle", "Mid-surface", "Surface")))

p17 <- ggplot(z_analysis, aes(x = z_layer, y = inst_speed_um_min, fill = z_layer)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.1, outlier.size = 0.3) +
  scale_fill_viridis_d() +
  coord_cartesian(ylim = c(0, quantile(z_analysis$inst_speed_um_min, 0.99, na.rm = TRUE))) +
  labs(title = "Nuclear Speed by Z-Depth",
       x = "Z-layer", y = "Instantaneous speed (µm/min)") +
  theme_track() +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "09_z_depth.pdf"), p17, width = 8, height = 5)
cat("  Saved: 09_z_depth.pdf\n")

# =============================================================================
# PART 12: PSEUDO-3D VISUALIZATION (ALL TRACKS)
# =============================================================================

cat("\n=== PART 12: PSEUDO-3D VISUALIZATION ===\n")

viz_all <- spots_clean %>%
  arrange(TRACK_ID, FRAME) %>%
  mutate(
    z_normalized    = (POSITION_Z - min(POSITION_Z)) / (max(POSITION_Z) - min(POSITION_Z) + 1e-10),
    time_normalized = FRAME / max(FRAME)
  )

# XY colored by Z
p_3d_xy <- ggplot(viz_all, aes(x = POSITION_X, y = POSITION_Y,
                                group = TRACK_ID, color = z_normalized)) +
  geom_path(alpha = 0.5, linewidth = 0.3) +
  scale_color_viridis_c(name = "Z depth\n(norm)", option = "plasma") +
  coord_fixed() +
  labs(title = "XY View (top-down)", subtitle = "Color = depth", x = "X (µm)", y = "Y (µm)") +
  theme_track()

# XY colored by time
p_3d_xy_time <- ggplot(viz_all, aes(x = POSITION_X, y = POSITION_Y,
                                     group = TRACK_ID, color = time_normalized)) +
  geom_path(alpha = 0.5, linewidth = 0.3) +
  scale_color_viridis_c(name = "Time\n(norm)", option = "viridis") +
  coord_fixed() +
  labs(title = "XY View", subtitle = "Color = time", x = "X (µm)", y = "Y (µm)") +
  theme_track()

# XZ side view
p_3d_xz <- ggplot(viz_all, aes(x = POSITION_X, y = POSITION_Z,
                                group = TRACK_ID, color = time_normalized)) +
  geom_path(alpha = 0.5, linewidth = 0.3) +
  scale_color_viridis_c(name = "Time", option = "viridis", guide = "none") +
  labs(title = "XZ View (side)", subtitle = "Color = time", x = "X (µm)", y = "Z (µm)") +
  theme_track()

# YZ front view
p_3d_yz <- ggplot(viz_all, aes(x = POSITION_Y, y = POSITION_Z,
                                group = TRACK_ID, color = time_normalized)) +
  geom_path(alpha = 0.5, linewidth = 0.3) +
  scale_color_viridis_c(name = "Time", option = "viridis", guide = "none") +
  labs(title = "YZ View (front)", subtitle = "Color = time", x = "Y (µm)", y = "Z (µm)") +
  theme_track()

pseudo_3d <- (p_3d_xy | p_3d_xy_time) / (p_3d_xz | p_3d_yz) +
  plot_annotation(
    title = "3D Track Visualization - All Tracks",
    subtitle = sprintf("%d tracks, %d spots | Multiple viewing angles",
                      n_distinct(viz_all$TRACK_ID), nrow(viz_all)),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey50")
    )
  )

ggsave(file.path(OUTPUT_DIR, "10_pseudo_3d_tracks.pdf"), pseudo_3d, width = 16, height = 14)
cat("  Saved: 10_pseudo_3d_tracks.pdf\n")

# =============================================================================
# PART 13: INTERACTIVE 3D/4D VISUALIZATION WITH PLOTLY
# =============================================================================

cat("\n=== PART 13: INTERACTIVE 4D VISUALIZATION ===\n")

# --- 13.1: 3D scatter by movement type ---
# Sample spots for interactive plotly (max 50k points for browser performance)
MAX_SPOTS_INTERACTIVE <- 50000
if (nrow(spots_clean) > MAX_SPOTS_INTERACTIVE) {
  set.seed(42)
  viz_interactive <- spots_clean %>%
    slice_sample(n = MAX_SPOTS_INTERACTIVE) %>%
    arrange(TRACK_ID, FRAME)
  cat(sprintf("  Sampled %d / %d spots for interactive plots\n",
              MAX_SPOTS_INTERACTIVE, nrow(spots_clean)))
} else {
  viz_interactive <- spots_clean %>% arrange(TRACK_ID, FRAME)
}

viz_interactive <- viz_interactive %>%
  left_join(
    tracks_classified %>% select(TRACK_ID, movement_type, SPEED_MEAN_UM_MIN),
    by = "TRACK_ID"
  ) %>%
  mutate(
    time_min = FRAME * FRAME_INTERVAL_MIN,
    movement_type = as.character(movement_type)
  )

fig_3d_movetype <- plot_ly(
  viz_interactive,
  x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
  color = ~movement_type,
  colors = c("Directed" = "#2166AC", "Random" = "#762A83", "Confined" = "#B2182B"),
  type = "scatter3d", mode = "markers",
  marker = list(size = 2, opacity = 0.5),
  hoverinfo = "text",
  text = ~paste("Track:", TRACK_ID, "<br>Type:", movement_type,
                "<br>Speed:", round(SPEED_MEAN_UM_MIN, 2), "µm/min")
) %>%
  plotly::layout(
    title = list(text = "3D Nuclear Positions by Movement Type", font = list(size = 16)),
    scene = list(
      xaxis = list(title = "X (µm)"),
      yaxis = list(title = "Y (µm)"),
      zaxis = list(title = "Z (µm)"),
      aspectmode = "data"
    )
  )

htmlwidgets::saveWidget(fig_3d_movetype,
                        file.path(OUTPUT_DIR, "11_interactive_3d_movetype.html"),
                        selfcontained = TRUE)
cat("  Saved: 11_interactive_3d_movetype.html\n")

# --- 13.2: 3D scatter by speed ---
fig_3d_speed <- plot_ly(
  viz_interactive,
  x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
  color = ~SPEED_MEAN_UM_MIN,
  colors = viridis::plasma(100),
  type = "scatter3d", mode = "markers",
  marker = list(size = 2, opacity = 0.5),
  hoverinfo = "text",
  text = ~paste("Track:", TRACK_ID, "<br>Speed:", round(SPEED_MEAN_UM_MIN, 2), "µm/min")
) %>%
  plotly::layout(
    title = list(text = "3D Nuclear Positions by Speed", font = list(size = 16)),
    scene = list(
      xaxis = list(title = "X (µm)"),
      yaxis = list(title = "Y (µm)"),
      zaxis = list(title = "Z (µm)"),
      aspectmode = "data"
    )
  )

htmlwidgets::saveWidget(fig_3d_speed,
                        file.path(OUTPUT_DIR, "12_interactive_3d_speed.html"),
                        selfcontained = TRUE)
cat("  Saved: 12_interactive_3d_speed.html\n")

# --- 13.3: 3D track lines (sample up to 500 tracks for performance) ---
MAX_TRACKS_3D <- 500
sampled_ids <- if (n_distinct(spots_clean$TRACK_ID) > MAX_TRACKS_3D) {
  sample(unique(spots_clean$TRACK_ID), MAX_TRACKS_3D)
} else {
  unique(spots_clean$TRACK_ID)
}
cat(sprintf("  Creating 3D track lines (%d tracks)...\n", length(sampled_ids)))

viz_lines <- spots_clean %>%
  filter(TRACK_ID %in% sampled_ids) %>%
  arrange(TRACK_ID, FRAME) %>%
  left_join(
    tracks_classified %>% select(TRACK_ID, movement_type, SPEED_MEAN_UM_MIN),
    by = "TRACK_ID"
  )

# Build traces grouped by movement type (avoids for-loop overhead)
fig_3d_lines <- plot_ly()

# Track which movement types have been added to legend
legend_shown <- c()

for (tid in unique(viz_lines$TRACK_ID)) {
  track_data <- viz_lines %>% filter(TRACK_ID == tid)
  mvtype <- as.character(track_data$movement_type[1])
  color <- switch(mvtype,
                  "Directed" = "#2166AC",
                  "Confined" = "#B2182B",
                  "#762A83")
  show_legend <- !(mvtype %in% legend_shown)
  if (show_legend) legend_shown <- c(legend_shown, mvtype)

  fig_3d_lines <- fig_3d_lines %>%
    add_trace(
      data = track_data,
      x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
      type = "scatter3d", mode = "lines",
      line = list(width = 3, color = color),
      opacity = 0.7,
      name = mvtype,
      legendgroup = mvtype,
      showlegend = show_legend,
      hoverinfo = "text",
      text = paste("Track:", tid, "| Type:", mvtype)
    )
}

fig_3d_lines <- fig_3d_lines %>%
  plotly::layout(
    title = list(text = sprintf("3D Track Trajectories (all %d tracks)",
                                n_distinct(viz_lines$TRACK_ID)),
                 font = list(size = 16)),
    scene = list(
      xaxis = list(title = "X (µm)"),
      yaxis = list(title = "Y (µm)"),
      zaxis = list(title = "Z (µm)"),
      aspectmode = "data"
    )
  )

htmlwidgets::saveWidget(fig_3d_lines,
                        file.path(OUTPUT_DIR, "13_interactive_3d_tracks.html"),
                        selfcontained = TRUE)
cat("  Saved: 13_interactive_3d_tracks.html\n")

# --- 13.4: 4D animation (time slider) ---
cat("  Creating 4D time-lapse animation...\n")

time_bins_anim <- seq(0, max(viz_interactive$time_min), by = 5)

viz_animated <- viz_interactive %>%
  mutate(time_bin = cut(time_min,
                        breaks = c(time_bins_anim, Inf),
                        labels = paste0(time_bins_anim, "-",
                                        c(time_bins_anim[-1], ceiling(max(time_min))), " min"),
                        include.lowest = TRUE))

fig_4d <- plot_ly(
  viz_animated,
  x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
  color = ~movement_type,
  colors = c("Directed" = "#2166AC", "Random" = "#762A83", "Confined" = "#B2182B"),
  frame = ~time_bin,
  type = "scatter3d", mode = "markers",
  marker = list(size = 2.5, opacity = 0.6)
) %>%
  plotly::layout(
    title = list(text = "4D Nuclear Movement (use slider for time)", font = list(size = 16)),
    scene = list(
      xaxis = list(title = "X (µm)"),
      yaxis = list(title = "Y (µm)"),
      zaxis = list(title = "Z (µm)"),
      aspectmode = "data"
    )
  ) %>%
  animation_opts(frame = 500, transition = 300, redraw = FALSE) %>%
  animation_slider(currentvalue = list(prefix = "Time: "))

htmlwidgets::saveWidget(fig_4d,
                        file.path(OUTPUT_DIR, "14_interactive_4d_timelapse.html"),
                        selfcontained = TRUE)
cat("  Saved: 14_interactive_4d_timelapse.html\n")

# =============================================================================
# PART 14: SPEED-CONFINEMENT PHASE SPACE
# =============================================================================

cat("\n=== PART 14: PHASE SPACE & SUMMARY ===\n")

p_phase <- ggplot(tracks_classified,
                  aes(x = SPEED_MEAN_UM_MIN, y = CONFINEMENT_RATIO, color = movement_type)) +
  geom_point(alpha = 0.6, size = 2.5) +
  scale_color_manual(values = movement_colors) +
  labs(title = "Speed-Confinement Phase Space",
       x = "Mean speed (µm/min)", y = "Confinement ratio (straightness)", color = NULL) +
  theme_track()

ggsave(file.path(OUTPUT_DIR, "15_phase_space.pdf"), p_phase, width = 8, height = 6)
cat("  Saved: 15_phase_space.pdf\n")

# =============================================================================
# EXPORT & FINAL SUMMARY
# =============================================================================

write_csv <- readr::write_csv
write_csv(tracks_classified, file.path(OUTPUT_DIR, "tracks_analyzed.csv"))
write_csv(spots_velocity,    file.path(OUTPUT_DIR, "spots_with_velocity.csv"))
cat("  Saved: tracks_analyzed.csv, spots_with_velocity.csv\n")

cat("\n")
cat(strrep("=", 70), "\n")
cat("IMARIS CSV ANALYSIS SUMMARY\n")
cat(strrep("=", 70), "\n\n")

cat(sprintf("DATASET: %s/ (prefix: %s)\n", DATA_DIR, FILE_PREFIX))
cat(sprintf("  Tracks: %d (kept %d after QC)\n", nrow(tracks), nrow(tracks_clean)))
cat(sprintf("  Spots:  %d (%d after QC)\n", nrow(spots), nrow(spots_clean)))
cat(sprintf("  Total imaging time: %.1f min (%.1f hours)\n",
            max(spots$FRAME) * FRAME_INTERVAL_MIN,
            max(spots$FRAME) * FRAME_INTERVAL_MIN / 60))

cat("\nNUCLEAR MOVEMENT:\n")
cat(sprintf("  Median speed: %.2f µm/min\n", median(tracks_clean$SPEED_MEAN_UM_MIN)))
cat(sprintf("  Median duration: %.1f min\n", median(tracks_clean$DURATION_MIN)))
cat(sprintf("  Median displacement: %.1f µm\n", median(tracks_clean$TRACK_DISPLACEMENT)))

cat("\nMOVEMENT TYPES:\n")
for (i in 1:nrow(movement_summary)) {
  cat(sprintf("  %s: %.1f%%\n", movement_summary$movement_type[i], movement_summary$pct[i]))
}

if (!is.na(alpha)) {
  cat(sprintf("\nMSD: α = %.2f → %s\n", alpha,
              ifelse(alpha < 0.8, "subdiffusive", ifelse(alpha > 1.2, "superdiffusive", "diffusive"))))
}

cat("\nOUTPUT FILES (in ", OUTPUT_DIR, "):\n")
cat("  01_quality_control.pdf\n")
cat("  02_velocity_analysis.pdf\n")
cat("  03_directionality.pdf\n")
cat("  04_spatial_velocity.pdf\n")
cat("  05_msd_analysis.pdf\n")
cat("  06_movement_classification.pdf\n")
cat("  08_confinement.pdf\n")
cat("  09_z_depth.pdf\n")
cat("  10_pseudo_3d_tracks.pdf\n")
cat("  11_interactive_3d_movetype.html\n")
cat("  12_interactive_3d_speed.html\n")
cat("  13_interactive_3d_tracks.html\n")
cat("  14_interactive_4d_timelapse.html\n")
cat("  15_phase_space.pdf\n")
cat("  tracks_analyzed.csv\n")
cat("  spots_with_velocity.csv\n")

cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n")
