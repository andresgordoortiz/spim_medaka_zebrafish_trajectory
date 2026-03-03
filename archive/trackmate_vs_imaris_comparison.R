# =============================================================================
# TrackMate vs Imaris — Quantitative & Qualitative Tracking Comparison
# =============================================================================
#
# PURPOSE:
#   Both TrackMate and Imaris were run on the SAME input images.
#   This script compares detection and tracking results to understand
#   how the two software packages differ in:
#     1. Spot detection (positions, counts per frame)
#     2. Track linking (number, duration, gaps)
#     3. Derived metrics (speed, displacement, confinement)
#     4. Movement classification
#     5. MSD / diffusion analysis
#
# DATA SOURCES:
#   TrackMate:  exported_from_trackmate_gui/4D_spots.csv & 4D_tracks.csv
#   Imaris:     all_tracks/ome-tiff.companion-bin-8bit-crop_*.csv
#
# CRITICAL NOTE ON TIME CALIBRATION:
#   Both software received images without proper time metadata.
#   The FRAME_INTERVAL must be set identically for a fair comparison.
#   If TrackMate and Imaris scripts used different intervals, speed/duration
#   values will be incomparable. Set FRAME_INTERVAL_SEC below to the TRUE
#   acquisition interval and both datasets will be re-calibrated here.
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
library(plotly)
library(htmlwidgets)
library(RANN)   # for nearest-neighbor spot matching

# =============================================================================
# GLOBAL PARAMETERS — SET THESE FOR YOUR EXPERIMENT
# =============================================================================

# *** USE A SINGLE FRAME INTERVAL FOR FAIR COMPARISON ***
# Both datasets come from the same images, so the true interval is the same.
# The original scripts used different values (TM: 30s, Imaris: 120s).
# Set this to your ACTUAL acquisition interval.
FRAME_INTERVAL_SEC <- 30
FRAME_INTERVAL_MIN <- FRAME_INTERVAL_SEC / 60

OUTPUT_DIR <- "analysis_output_comparison"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# QC thresholds (applied identically to both datasets)
MIN_SPOTS <- 5
MAX_SPEED_UM_MIN <- 10

# =============================================================================
# PLOT STYLING
# =============================================================================

software_colors <- c("TrackMate" = "#2166AC", "Imaris" = "#B2182B")

movement_colors <- c(
  "Directed" = "#2166AC",
  "Confined" = "#B2182B",
  "Random"   = "#762A83"
)

theme_compare <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1, color = "grey50"),
      legend.position = "bottom",
      legend.key.size = unit(0.3, "cm"),
      legend.text     = element_text(size = base_size - 1),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = base_size - 1),
      strip.text = element_text(face = "bold", size = base_size)
    )
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("  TRACKMATE vs IMARIS — TRACKING COMPARISON\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("  Frame interval for comparison: %d sec (%.2f min)\n",
            FRAME_INTERVAL_SEC, FRAME_INTERVAL_MIN))

# #############################################################################
# PART 1: DATA LOADING
# #############################################################################

cat("\n=== PART 1: DATA LOADING ===\n")

# --- 1A: TrackMate ---
cat("  Loading TrackMate data...\n")

tm_spots_raw <- read_csv("exported_from_trackmate_gui/4D_spots.csv", show_col_types = FALSE)[-c(1:3), ]
tm_tracks_raw <- read_csv("exported_from_trackmate_gui/4D_tracks.csv", show_col_types = FALSE)[-c(1:3), ]

tm_spots <- tm_spots_raw %>%
  mutate(across(c(ID, TRACK_ID, QUALITY, POSITION_X, POSITION_Y, POSITION_Z,
                  POSITION_T, FRAME, RADIUS),
                ~as.numeric(.))) %>%
  filter(!is.na(TRACK_ID)) %>%
  select(ID, TRACK_ID, POSITION_X, POSITION_Y, POSITION_Z, FRAME, QUALITY)

tm_tracks <- tm_tracks_raw %>%
  mutate(across(-LABEL, ~as.numeric(.))) %>%
  mutate(
    DURATION_MIN      = TRACK_DURATION * FRAME_INTERVAL_MIN,
    SPEED_MEAN_UM_MIN = TRACK_MEAN_SPEED * (60 / FRAME_INTERVAL_SEC),
    N_SPOTS           = NUMBER_SPOTS
  )

# --- 1B: Imaris ---
cat("  Loading Imaris data...\n")

DATA_DIR_IM   <- "all_tracks"
FILE_PREFIX   <- "ome-tiff.companion-bin-8bit-crop"

read_imaris_csv <- function(stat_name) {
  path <- file.path(DATA_DIR_IM, paste0(FILE_PREFIX, "_", stat_name, ".csv"))
  df <- read_csv(path, skip = 3, show_col_types = FALSE)
  df %>% select(where(~ !all(is.na(.x))))
}

im_spots <- read_imaris_csv("Position") %>%
  rename(POSITION_X = `Position X`, POSITION_Y = `Position Y`,
         POSITION_Z = `Position Z`, FRAME = Time,
         TRACK_ID = TrackID, ID = ID) %>%
  select(ID, TRACK_ID, POSITION_X, POSITION_Y, POSITION_Z, FRAME) %>%
  mutate(across(everything(), as.numeric))

im_track_dur <- read_imaris_csv("Track_Duration") %>%
  rename(TRACK_DURATION_FRAMES = `Track Duration`, TRACK_ID = ID) %>%
  select(TRACK_ID, TRACK_DURATION_FRAMES) %>%
  mutate(across(everything(), as.numeric))

im_track_len <- read_imaris_csv("Track_Length") %>%
  rename(TRACK_LENGTH = `Track Length`, TRACK_ID = ID) %>%
  select(TRACK_ID, TRACK_LENGTH) %>%
  mutate(across(everything(), as.numeric))

im_track_spd <- read_imaris_csv("Track_Speed_Mean") %>%
  rename(TRACK_SPEED_MEAN_RAW = `Track Speed Mean`, TRACK_ID = ID) %>%
  select(TRACK_ID, TRACK_SPEED_MEAN_RAW) %>%
  mutate(across(everything(), as.numeric))

im_track_str <- read_imaris_csv("Track_Straightness") %>%
  rename(CONFINEMENT_RATIO = `Track Straightness`, TRACK_ID = ID) %>%
  select(TRACK_ID, CONFINEMENT_RATIO) %>%
  mutate(across(everything(), as.numeric))

im_track_disp <- read_imaris_csv("Track_Displacement_Length") %>%
  rename(TRACK_DISPLACEMENT = `Track Displacement Length`, TRACK_ID = ID) %>%
  select(TRACK_ID, TRACK_DISPLACEMENT) %>%
  mutate(across(everything(), as.numeric))

im_spots_per_track <- im_spots %>% count(TRACK_ID, name = "N_SPOTS")

im_tracks <- im_track_dur %>%
  left_join(im_track_len,  by = "TRACK_ID") %>%
  left_join(im_track_spd,  by = "TRACK_ID") %>%
  left_join(im_track_str,  by = "TRACK_ID") %>%
  left_join(im_track_disp, by = "TRACK_ID") %>%
  left_join(im_spots_per_track, by = "TRACK_ID") %>%
  mutate(
    # Imaris "duration" = number of frames (uncalibrated as 1 s/frame)
    DURATION_MIN      = TRACK_DURATION_FRAMES * FRAME_INTERVAL_SEC / 60,
    # Imaris "speed" = µm/frame (uncalibrated time)
    SPEED_MEAN_UM_MIN = TRACK_SPEED_MEAN_RAW / FRAME_INTERVAL_SEC * 60
  )

# --- Print raw counts ---
cat(sprintf("\n  %-12s  %10s  %10s\n", "", "TrackMate", "Imaris"))
cat(sprintf("  %-12s  %10d  %10d\n", "Spots",  nrow(tm_spots),  nrow(im_spots)))
cat(sprintf("  %-12s  %10d  %10d\n", "Tracks", nrow(tm_tracks), nrow(im_tracks)))
cat(sprintf("  %-12s  %10d  %10d\n", "Max frame",
            max(tm_spots$FRAME, na.rm = TRUE), max(im_spots$FRAME, na.rm = TRUE)))

# #############################################################################
# PART 2: QUALITY CONTROL (same criteria for both)
# #############################################################################

cat("\n=== PART 2: QUALITY CONTROL (uniform criteria) ===\n")
cat(sprintf("  MIN_SPOTS = %d | MAX_SPEED = %d µm/min\n", MIN_SPOTS, MAX_SPEED_UM_MIN))

tm_tracks <- tm_tracks %>%
  mutate(passes_qc = N_SPOTS >= MIN_SPOTS & SPEED_MEAN_UM_MIN <= MAX_SPEED_UM_MIN)
im_tracks <- im_tracks %>%
  mutate(passes_qc = N_SPOTS >= MIN_SPOTS & SPEED_MEAN_UM_MIN <= MAX_SPEED_UM_MIN)

tm_tracks_clean <- tm_tracks %>% filter(passes_qc)
im_tracks_clean <- im_tracks %>% filter(passes_qc)

tm_spots_clean <- tm_spots %>% filter(TRACK_ID %in% tm_tracks_clean$TRACK_ID)
im_spots_clean <- im_spots %>% filter(TRACK_ID %in% im_tracks_clean$TRACK_ID)

cat(sprintf("\n  After QC:    TrackMate    Imaris\n"))
cat(sprintf("  Tracks:      %7d     %7d\n", nrow(tm_tracks_clean), nrow(im_tracks_clean)))
cat(sprintf("  Spots:       %7d     %7d\n", nrow(tm_spots_clean),  nrow(im_spots_clean)))

# #############################################################################
# PART 3: DETECTION COMPARISON
# #############################################################################

cat("\n=== PART 3: DETECTION COMPARISON ===\n")

# --- 3A: Spots per frame ---
tm_spf <- tm_spots %>%
  count(FRAME, name = "n_spots") %>%
  mutate(software = "TrackMate")

im_spf <- im_spots %>%
  count(FRAME, name = "n_spots") %>%
  mutate(software = "Imaris")

spots_per_frame <- bind_rows(tm_spf, im_spf) %>%
  mutate(time_min = FRAME * FRAME_INTERVAL_MIN)

p_spf <- ggplot(spots_per_frame, aes(x = time_min, y = n_spots, color = software)) +
  geom_line(alpha = 0.8, linewidth = 0.5) +
  scale_color_manual(values = software_colors) +
  labs(title = "Spots Detected per Frame",
       subtitle = "Same images, different detection algorithms",
       x = "Time (min)", y = "Number of spots", color = NULL) +
  theme_compare()

# --- 3B: Spot position distributions (marginal) ---
tm_pos <- tm_spots %>% mutate(software = "TrackMate")
im_pos <- im_spots %>% mutate(software = "Imaris")
all_pos <- bind_rows(
  tm_pos %>% select(POSITION_X, POSITION_Y, POSITION_Z, software),
  im_pos %>% select(POSITION_X, POSITION_Y, POSITION_Z, software)
)

p_xdist <- ggplot(all_pos, aes(x = POSITION_X, fill = software)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = software_colors) +
  labs(title = "X Position Distribution", x = "X (µm)", y = "Density", fill = NULL) +
  theme_compare()

p_ydist <- ggplot(all_pos, aes(x = POSITION_Y, fill = software)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = software_colors) +
  labs(title = "Y Position Distribution", x = "Y (µm)", y = "Density", fill = NULL) +
  theme_compare()

p_zdist <- ggplot(all_pos, aes(x = POSITION_Z, fill = software)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = software_colors) +
  labs(title = "Z Position Distribution", x = "Z (µm)", y = "Density", fill = NULL) +
  theme_compare()

detection_combined <- p_spf / (p_xdist + p_ydist + p_zdist) +
  plot_annotation(
    title = "Detection Comparison: TrackMate vs Imaris",
    subtitle = "Same input images — comparing spot detection",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey50"))
  )

ggsave(file.path(OUTPUT_DIR, "01_detection_comparison.pdf"), detection_combined,
       width = 14, height = 10)
cat("  Saved: 01_detection_comparison.pdf\n")

# #############################################################################
# PART 4: SPOT-LEVEL SPATIAL MATCHING (per-frame nearest neighbor)
# #############################################################################

cat("\n=== PART 4: SPOT MATCHING (nearest-neighbor) ===\n")
cat("  Matching spots frame-by-frame using 3D nearest neighbor...\n")

# Sample frames for spot matching (doing all frames is expensive)
shared_frames <- intersect(unique(tm_spots$FRAME), unique(im_spots$FRAME))
set.seed(42)
sample_frames <- sort(sample(shared_frames, min(50, length(shared_frames))))

match_results <- map_dfr(sample_frames, function(f) {
  tm_f <- tm_spots %>% filter(FRAME == f) %>% select(POSITION_X, POSITION_Y, POSITION_Z)
  im_f <- im_spots %>% filter(FRAME == f) %>% select(POSITION_X, POSITION_Y, POSITION_Z)

  if (nrow(tm_f) < 2 || nrow(im_f) < 2) return(NULL)

  # For each TrackMate spot, find nearest Imaris spot
  nn <- RANN::nn2(as.matrix(im_f), as.matrix(tm_f), k = 1)
  tibble(
    FRAME = f,
    nn_dist = as.numeric(nn$nn.dists),
    direction = "TM→IM"
  )
})

# Also in reverse: for each Imaris spot, find nearest TrackMate spot
match_results_rev <- map_dfr(sample_frames, function(f) {
  tm_f <- tm_spots %>% filter(FRAME == f) %>% select(POSITION_X, POSITION_Y, POSITION_Z)
  im_f <- im_spots %>% filter(FRAME == f) %>% select(POSITION_X, POSITION_Y, POSITION_Z)

  if (nrow(tm_f) < 2 || nrow(im_f) < 2) return(NULL)

  nn <- RANN::nn2(as.matrix(tm_f), as.matrix(im_f), k = 1)
  tibble(
    FRAME = f,
    nn_dist = as.numeric(nn$nn.dists),
    direction = "IM→TM"
  )
})

match_all <- bind_rows(match_results, match_results_rev)

# Thresholds for "matched"
MATCH_THRESHOLD_UM <- 5  # µm — spots within 5 µm are considered the same nucleus

match_summary <- match_all %>%
  group_by(direction) %>%
  summarise(
    total          = n(),
    matched        = sum(nn_dist < MATCH_THRESHOLD_UM),
    pct_matched    = 100 * matched / total,
    median_dist    = median(nn_dist),
    mean_dist      = mean(nn_dist),
    q95_dist       = quantile(nn_dist, 0.95),
    .groups = "drop"
  )

cat("\n  Nearest-Neighbor Matching (sampled frames):\n")
cat(sprintf("  %-7s  matched  pct_matched  median_dist  mean_dist  q95_dist\n", "dir"))
for (i in 1:nrow(match_summary)) {
  cat(sprintf("  %-7s  %7d  %10.1f%%  %10.2f µm  %8.2f µm  %7.2f µm\n",
              match_summary$direction[i], match_summary$matched[i],
              match_summary$pct_matched[i], match_summary$median_dist[i],
              match_summary$mean_dist[i], match_summary$q95_dist[i]))
}

# Plot: distance to nearest match
p_match <- ggplot(match_all, aes(x = nn_dist, fill = direction)) +
  geom_histogram(bins = 100, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = MATCH_THRESHOLD_UM, linetype = "dashed", color = "grey30") +
  scale_fill_manual(values = c("TM→IM" = "#2166AC", "IM→TM" = "#B2182B")) +
  scale_x_continuous(limits = c(0, 50)) +
  labs(title = "Distance to Nearest Matched Spot",
       subtitle = sprintf("Dashed = %d µm threshold | Per-frame 3D nearest neighbor", MATCH_THRESHOLD_UM),
       x = "Distance (µm)", y = "Count", fill = "Direction") +
  theme_compare()

# Matched fraction over time
match_over_time <- match_all %>%
  mutate(time_min = FRAME * FRAME_INTERVAL_MIN) %>%
  group_by(time_min, direction) %>%
  summarise(pct_matched = 100 * mean(nn_dist < MATCH_THRESHOLD_UM), .groups = "drop")

p_match_time <- ggplot(match_over_time, aes(x = time_min, y = pct_matched, color = direction)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("TM→IM" = "#2166AC", "IM→TM" = "#B2182B")) +
  labs(title = "Spot Matching Rate Over Time",
       subtitle = sprintf("Fraction of spots with a match within %d µm", MATCH_THRESHOLD_UM),
       x = "Time (min)", y = "% matched", color = "Direction") +
  theme_compare()

matching_combined <- (p_match + p_match_time) +
  plot_annotation(
    title = "Spatial Spot Matching: TrackMate vs Imaris",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "02_spot_matching.pdf"), matching_combined, width = 14, height = 5)
cat("  Saved: 02_spot_matching.pdf\n")

# #############################################################################
# PART 5: TRACK-LEVEL COMPARISON (distributions)
# #############################################################################

cat("\n=== PART 5: TRACK-LEVEL COMPARISON ===\n")

# Combine track data
tm_track_df <- tm_tracks_clean %>%
  transmute(
    software           = "TrackMate",
    N_SPOTS            = N_SPOTS,
    DURATION_MIN       = DURATION_MIN,
    SPEED_MEAN_UM_MIN  = SPEED_MEAN_UM_MIN,
    TRACK_DISPLACEMENT = TRACK_DISPLACEMENT,
    CONFINEMENT_RATIO  = CONFINEMENT_RATIO,
    TRACK_LENGTH       = TOTAL_DISTANCE_TRAVELED
  )

im_track_df <- im_tracks_clean %>%
  transmute(
    software           = "Imaris",
    N_SPOTS            = N_SPOTS,
    DURATION_MIN       = DURATION_MIN,
    SPEED_MEAN_UM_MIN  = SPEED_MEAN_UM_MIN,
    TRACK_DISPLACEMENT = TRACK_DISPLACEMENT,
    CONFINEMENT_RATIO  = CONFINEMENT_RATIO,
    TRACK_LENGTH       = TRACK_LENGTH
  )

all_tracks <- bind_rows(tm_track_df, im_track_df) %>%
  mutate(software = factor(software, levels = c("TrackMate", "Imaris")))

# --- 5A: Number of spots per track ---
p_nspots <- ggplot(all_tracks, aes(x = N_SPOTS, fill = software)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = software_colors) +
  scale_x_log10() +
  labs(title = "Track Length (spots per track)", x = "N spots (log)", y = "Count", fill = NULL) +
  theme_compare()

# --- 5B: Duration ---
p_dur <- ggplot(all_tracks, aes(x = DURATION_MIN, fill = software)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = software_colors) +
  scale_x_log10() +
  labs(title = "Track Duration", x = "Duration (min, log)", y = "Count", fill = NULL) +
  theme_compare()

# --- 5C: Mean speed ---
p_speed <- ggplot(all_tracks, aes(x = SPEED_MEAN_UM_MIN, fill = software)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = software_colors) +
  coord_cartesian(xlim = c(0, quantile(all_tracks$SPEED_MEAN_UM_MIN, 0.99, na.rm = TRUE))) +
  labs(title = "Mean Track Speed", x = "Speed (µm/min)", y = "Density", fill = NULL) +
  theme_compare()

# --- 5D: Displacement ---
p_disp <- ggplot(all_tracks, aes(x = TRACK_DISPLACEMENT, fill = software)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = software_colors) +
  coord_cartesian(xlim = c(0, quantile(all_tracks$TRACK_DISPLACEMENT, 0.99, na.rm = TRUE))) +
  labs(title = "Track Displacement", x = "Displacement (µm)", y = "Density", fill = NULL) +
  theme_compare()

# --- 5E: Confinement ratio ---
p_confine <- ggplot(all_tracks, aes(x = CONFINEMENT_RATIO, fill = software)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = software_colors) +
  labs(title = "Confinement Ratio (Straightness)",
       x = "Confinement ratio", y = "Density", fill = NULL) +
  theme_compare()

# --- 5F: Track length (total distance traveled) ---
p_pathlen <- ggplot(all_tracks, aes(x = TRACK_LENGTH, fill = software)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = software_colors) +
  coord_cartesian(xlim = c(0, quantile(all_tracks$TRACK_LENGTH, 0.99, na.rm = TRUE))) +
  labs(title = "Total Path Length", x = "Path length (µm)", y = "Density", fill = NULL) +
  theme_compare()

track_combined <- (p_nspots + p_dur) / (p_speed + p_disp) / (p_confine + p_pathlen) +
  plot_annotation(
    title = "Track-Level Metric Comparison: TrackMate vs Imaris",
    subtitle = sprintf("TrackMate: %d tracks | Imaris: %d tracks (after QC)",
                       nrow(tm_tracks_clean), nrow(im_tracks_clean)),
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey50"))
  )

ggsave(file.path(OUTPUT_DIR, "03_track_distributions.pdf"), track_combined,
       width = 14, height = 14)
cat("  Saved: 03_track_distributions.pdf\n")

# --- Summary statistics table ---
track_stats_table <- all_tracks %>%
  group_by(software) %>%
  summarise(
    n_tracks             = n(),
    median_n_spots       = median(N_SPOTS, na.rm = TRUE),
    median_duration_min  = median(DURATION_MIN, na.rm = TRUE),
    median_speed_um_min  = median(SPEED_MEAN_UM_MIN, na.rm = TRUE),
    mean_speed_um_min    = mean(SPEED_MEAN_UM_MIN, na.rm = TRUE),
    median_displacement  = median(TRACK_DISPLACEMENT, na.rm = TRUE),
    median_confinement   = median(CONFINEMENT_RATIO, na.rm = TRUE),
    median_path_length   = median(TRACK_LENGTH, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n  Track Summary Statistics:\n")
print(track_stats_table, width = Inf)

# --- Statistical tests ---
cat("\n  Statistical Tests (Wilcoxon rank-sum):\n")
test_vars <- c("N_SPOTS", "DURATION_MIN", "SPEED_MEAN_UM_MIN",
               "TRACK_DISPLACEMENT", "CONFINEMENT_RATIO", "TRACK_LENGTH")

for (v in test_vars) {
  tm_vals <- tm_track_df[[v]]
  im_vals <- im_track_df[[v]]
  wt <- wilcox.test(tm_vals, im_vals, conf.int = FALSE)
  effect_size <- abs(median(tm_vals, na.rm = TRUE) - median(im_vals, na.rm = TRUE)) /
    max(median(tm_vals, na.rm = TRUE), median(im_vals, na.rm = TRUE))
  cat(sprintf("  %-22s  p = %.2e  |  rel. diff = %.1f%%\n",
              v, wt$p.value, 100 * effect_size))
}

# #############################################################################
# PART 6: VELOCITY ANALYSIS COMPARISON
# #############################################################################

cat("\n=== PART 6: VELOCITY COMPARISON ===\n")

compute_velocities <- function(spots_df, frame_int_min) {
  spots_df %>%
    arrange(TRACK_ID, FRAME) %>%
    group_by(TRACK_ID) %>%
    mutate(
      dx = POSITION_X - lag(POSITION_X),
      dy = POSITION_Y - lag(POSITION_Y),
      dz = POSITION_Z - lag(POSITION_Z),
      dt_min = (FRAME - lag(FRAME)) * frame_int_min,
      displacement = sqrt(dx^2 + dy^2 + dz^2),
      inst_speed = displacement / dt_min,
      vx = dx / dt_min, vy = dy / dt_min, vz = dz / dt_min,
      time_min = FRAME * frame_int_min
    ) %>%
    ungroup() %>%
    filter(!is.na(inst_speed) & dt_min > 0)
}

tm_vel <- compute_velocities(tm_spots_clean, FRAME_INTERVAL_MIN) %>%
  mutate(software = "TrackMate")
im_vel <- compute_velocities(im_spots_clean, FRAME_INTERVAL_MIN) %>%
  mutate(software = "Imaris")

all_vel <- bind_rows(tm_vel, im_vel)

# Instantaneous speed distributions
p_inst_speed <- ggplot(all_vel, aes(x = inst_speed, fill = software)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = software_colors) +
  coord_cartesian(xlim = c(0, quantile(all_vel$inst_speed, 0.99, na.rm = TRUE))) +
  labs(title = "Instantaneous Speed Distribution",
       x = "Speed (µm/min)", y = "Density", fill = NULL) +
  theme_compare()

# Speed over time (binned)
speed_time <- all_vel %>%
  mutate(time_bin = floor(time_min / 10) * 10) %>%
  group_by(time_bin, software) %>%
  summarise(
    mean_speed = mean(inst_speed, na.rm = TRUE),
    se_speed   = sd(inst_speed, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_speed_time <- ggplot(speed_time, aes(x = time_bin, y = mean_speed, color = software)) +
  geom_ribbon(aes(ymin = mean_speed - se_speed, ymax = mean_speed + se_speed,
                  fill = software), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = software_colors) +
  scale_fill_manual(values = software_colors) +
  labs(title = "Mean Instantaneous Speed Over Time",
       x = "Time (min)", y = "Mean speed (µm/min)", color = NULL, fill = NULL) +
  theme_compare()

# Turning angle comparison
compute_turning <- function(vel_df) {
  vel_df %>%
    filter(!is.na(dx) & !is.na(dy) & !is.na(dz)) %>%
    group_by(TRACK_ID) %>%
    mutate(
      dot = dx * lag(dx) + dy * lag(dy) + dz * lag(dz),
      mag1 = sqrt(lag(dx)^2 + lag(dy)^2 + lag(dz)^2),
      mag2 = sqrt(dx^2 + dy^2 + dz^2),
      turning_angle = acos(pmin(pmax(dot / (mag1 * mag2 + 1e-12), -1), 1)) * 180 / pi
    ) %>%
    ungroup() %>%
    filter(!is.na(turning_angle))
}

tm_turn <- compute_turning(tm_vel) %>% mutate(software = "TrackMate")
im_turn <- compute_turning(im_vel) %>% mutate(software = "Imaris")
all_turn <- bind_rows(tm_turn, im_turn)

p_turning <- ggplot(all_turn, aes(x = turning_angle, fill = software)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = software_colors) +
  geom_vline(xintercept = 90, linetype = "dashed", color = "grey50") +
  labs(title = "Turning Angle Distribution",
       subtitle = "90° = random; <90° = persistent; >90° = reversing",
       x = "Turning angle (°)", y = "Density", fill = NULL) +
  theme_compare()

velocity_combined <- (p_inst_speed + p_speed_time) / (p_turning + plot_spacer()) +
  plot_annotation(
    title = "Velocity & Directionality Comparison",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "04_velocity_comparison.pdf"), velocity_combined,
       width = 14, height = 10)
cat("  Saved: 04_velocity_comparison.pdf\n")

# #############################################################################
# PART 7: MSD COMPARISON
# #############################################################################

cat("\n=== PART 7: MSD COMPARISON ===\n")

calculate_msd <- function(df, max_lag = 30) {
  df <- df %>% arrange(FRAME)
  n <- nrow(df)
  if (n < 5) return(NULL)
  max_lag <- min(max_lag, n - 1)
  msd_values <- sapply(1:max_lag, function(lag) {
    d <- (df$POSITION_X[(lag+1):n] - df$POSITION_X[1:(n-lag)])^2 +
         (df$POSITION_Y[(lag+1):n] - df$POSITION_Y[1:(n-lag)])^2 +
         (df$POSITION_Z[(lag+1):n] - df$POSITION_Z[1:(n-lag)])^2
    mean(d, na.rm = TRUE)
  })
  tibble(lag_frames = 1:max_lag, lag_min = 1:max_lag * FRAME_INTERVAL_MIN, msd = msd_values)
}

MAX_MSD_TRACKS <- 500

compute_ensemble_msd <- function(spots_df, sw_name) {
  track_ids <- unique(spots_df$TRACK_ID)
  if (length(track_ids) > MAX_MSD_TRACKS) {
    set.seed(42)
    track_ids <- sample(track_ids, MAX_MSD_TRACKS)
  }
  cat(sprintf("  Computing MSD for %s (%d tracks)...\n", sw_name, length(track_ids)))

  msd_data <- spots_df %>%
    filter(TRACK_ID %in% track_ids) %>%
    arrange(TRACK_ID, FRAME) %>%
    group_by(TRACK_ID) %>%
    group_modify(~{
      result <- calculate_msd(.x)
      if (is.null(result)) return(tibble())
      result
    }) %>%
    ungroup()

  ensemble <- msd_data %>%
    group_by(lag_min) %>%
    summarise(
      mean_msd = mean(msd), se_msd = sd(msd) / sqrt(n()), n = n(),
      .groups = "drop"
    ) %>%
    filter(n >= 10) %>%
    mutate(software = sw_name)

  # Fit alpha
  fit <- tryCatch(
    lm(log(mean_msd) ~ log(lag_min), data = ensemble %>% filter(lag_min <= 10)),
    error = function(e) NULL
  )
  alpha <- if (!is.null(fit)) coef(fit)[2] else NA
  cat(sprintf("  %s α = %.3f\n", sw_name, alpha))

  list(individual = msd_data %>% mutate(software = sw_name),
       ensemble = ensemble,
       alpha = alpha)
}

tm_msd <- compute_ensemble_msd(tm_spots_clean, "TrackMate")
im_msd <- compute_ensemble_msd(im_spots_clean, "Imaris")

ensemble_both <- bind_rows(tm_msd$ensemble, im_msd$ensemble)

p_msd <- ggplot(ensemble_both, aes(x = lag_min, y = mean_msd, color = software)) +
  geom_ribbon(aes(ymin = mean_msd - se_msd, ymax = mean_msd + se_msd, fill = software),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = software_colors) +
  scale_fill_manual(values = software_colors) +
  scale_x_log10() + scale_y_log10() +
  labs(
    title = "Mean Squared Displacement Comparison",
    subtitle = sprintf("TrackMate α = %.2f | Imaris α = %.2f",
                       tm_msd$alpha, im_msd$alpha),
    x = "Time lag (min, log)", y = "MSD (µm², log)", color = NULL, fill = NULL
  ) +
  theme_compare()

ggsave(file.path(OUTPUT_DIR, "05_msd_comparison.pdf"), p_msd, width = 8, height = 6)
cat("  Saved: 05_msd_comparison.pdf\n")

# #############################################################################
# PART 8: MOVEMENT TYPE CLASSIFICATION COMPARISON
# #############################################################################

cat("\n=== PART 8: MOVEMENT CLASSIFICATION COMPARISON ===\n")

# Use a uniform classification for both (confinement + speed based)
classify_tracks <- function(track_df, sw_name) {
  track_df %>%
    mutate(
      movement_type = case_when(
        CONFINEMENT_RATIO > 0.5 & SPEED_MEAN_UM_MIN > median(SPEED_MEAN_UM_MIN, na.rm = TRUE) ~ "Directed",
        CONFINEMENT_RATIO < 0.25 ~ "Confined",
        TRUE ~ "Random"
      ),
      movement_type = factor(movement_type, levels = c("Directed", "Random", "Confined")),
      software = sw_name
    )
}

tm_classified <- classify_tracks(tm_tracks_clean, "TrackMate")
im_classified <- classify_tracks(im_tracks_clean, "Imaris")
all_classified <- bind_rows(
  tm_classified %>% select(software, movement_type, CONFINEMENT_RATIO, SPEED_MEAN_UM_MIN),
  im_classified %>% select(software, movement_type, CONFINEMENT_RATIO, SPEED_MEAN_UM_MIN)
)

move_summary <- all_classified %>%
  group_by(software, movement_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(software) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

cat("\n  Movement Type Distribution:\n")
cat(sprintf("  %-10s  %-10s  %6s  %6s\n", "Software", "Type", "N", "%"))
for (i in 1:nrow(move_summary)) {
  cat(sprintf("  %-10s  %-10s  %6d  %5.1f%%\n",
              move_summary$software[i], move_summary$movement_type[i],
              move_summary$n[i], move_summary$pct[i]))
}

# Grouped bar chart
p_move_bar <- ggplot(move_summary, aes(x = movement_type, y = pct, fill = software)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 3) +
  scale_fill_manual(values = software_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Movement Type Proportions",
       x = NULL, y = "% of tracks", fill = NULL) +
  theme_compare() +
  theme(panel.grid.major.x = element_blank())

# Phase space: confinement vs speed
p_phase <- ggplot(all_classified,
                  aes(x = CONFINEMENT_RATIO, y = SPEED_MEAN_UM_MIN, color = software)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_density2d(linewidth = 0.5) +
  scale_color_manual(values = software_colors) +
  facet_wrap(~software) +
  labs(title = "Speed-Confinement Phase Space",
       x = "Confinement ratio", y = "Speed (µm/min)", color = NULL) +
  theme_compare() +
  theme(legend.position = "none")

classification_combined <- (p_move_bar + p_phase) +
  plot_annotation(
    title = "Movement Classification Comparison",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "06_movement_classification.pdf"), classification_combined,
       width = 14, height = 6)
cat("  Saved: 06_movement_classification.pdf\n")

# #############################################################################
# PART 9: SPATIAL COMPARISON — SIDE-BY-SIDE TRACK PROJECTIONS
# #############################################################################

cat("\n=== PART 9: SPATIAL COMPARISON ===\n")

# Sample tracks for visualization
MAX_VIZ_TRACKS <- 500

sample_tracks <- function(spots_df, n_tracks) {
  ids <- unique(spots_df$TRACK_ID)
  if (length(ids) > n_tracks) {
    set.seed(42)
    ids <- sample(ids, n_tracks)
  }
  spots_df %>% filter(TRACK_ID %in% ids) %>% arrange(TRACK_ID, FRAME)
}

tm_viz <- sample_tracks(tm_spots_clean, MAX_VIZ_TRACKS) %>%
  mutate(z_norm = (POSITION_Z - min(POSITION_Z)) / (max(POSITION_Z) - min(POSITION_Z) + 1e-10))
im_viz <- sample_tracks(im_spots_clean, MAX_VIZ_TRACKS) %>%
  mutate(z_norm = (POSITION_Z - min(POSITION_Z)) / (max(POSITION_Z) - min(POSITION_Z) + 1e-10))

p_xy_tm <- ggplot(tm_viz, aes(x = POSITION_X, y = POSITION_Y,
                               group = TRACK_ID, color = z_norm)) +
  geom_path(alpha = 0.4, linewidth = 0.3) +
  scale_color_viridis_c(name = "Z depth", option = "plasma") +
  coord_fixed() +
  labs(title = "TrackMate — XY Projection", x = "X (µm)", y = "Y (µm)") +
  theme_compare()

p_xy_im <- ggplot(im_viz, aes(x = POSITION_X, y = POSITION_Y,
                               group = TRACK_ID, color = z_norm)) +
  geom_path(alpha = 0.4, linewidth = 0.3) +
  scale_color_viridis_c(name = "Z depth", option = "plasma") +
  coord_fixed() +
  labs(title = "Imaris — XY Projection", x = "X (µm)", y = "Y (µm)") +
  theme_compare()

p_xz_tm <- ggplot(tm_viz, aes(x = POSITION_X, y = POSITION_Z,
                               group = TRACK_ID, color = z_norm)) +
  geom_path(alpha = 0.4, linewidth = 0.3) +
  scale_color_viridis_c(guide = "none", option = "plasma") +
  labs(title = "TrackMate — XZ", x = "X (µm)", y = "Z (µm)") +
  theme_compare()

p_xz_im <- ggplot(im_viz, aes(x = POSITION_X, y = POSITION_Z,
                               group = TRACK_ID, color = z_norm)) +
  geom_path(alpha = 0.4, linewidth = 0.3) +
  scale_color_viridis_c(guide = "none", option = "plasma") +
  labs(title = "Imaris — XZ", x = "X (µm)", y = "Z (µm)") +
  theme_compare()

spatial_combined <- (p_xy_tm + p_xy_im) / (p_xz_tm + p_xz_im) +
  plot_annotation(
    title = "Spatial Track Comparison (500 sampled tracks each)",
    subtitle = "Same input images — different tracking algorithms",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey50"))
  )

ggsave(file.path(OUTPUT_DIR, "07_spatial_tracks.pdf"), spatial_combined,
       width = 16, height = 14)
cat("  Saved: 07_spatial_tracks.pdf\n")

# #############################################################################
# PART 10: Z-DEPTH COMPARISON
# #############################################################################

cat("\n=== PART 10: Z-DEPTH COMPARISON ===\n")

z_vel <- all_vel %>%
  mutate(z_layer = cut(POSITION_Z, breaks = 5,
                       labels = c("Deep", "Mid-deep", "Middle", "Mid-surface", "Surface")))

p_z_speed <- ggplot(z_vel, aes(x = z_layer, y = inst_speed, fill = software)) +
  geom_boxplot(alpha = 0.6, outlier.alpha = 0.05, outlier.size = 0.2,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = software_colors) +
  coord_cartesian(ylim = c(0, quantile(z_vel$inst_speed, 0.95, na.rm = TRUE))) +
  labs(title = "Instantaneous Speed by Z-Depth",
       x = "Z-layer", y = "Speed (µm/min)", fill = NULL) +
  theme_compare()

ggsave(file.path(OUTPUT_DIR, "08_z_depth_comparison.pdf"), p_z_speed, width = 10, height = 6)
cat("  Saved: 08_z_depth_comparison.pdf\n")

# #############################################################################
# PART 11: SPATIAL SPEED HEATMAPS COMPARISON
# #############################################################################

cat("\n=== PART 11: SPATIAL SPEED HEATMAPS ===\n")

BIN_SIZE <- 30  # µm

make_heatmap_data <- function(vel_df, bin_size) {
  vel_df %>%
    mutate(x_bin = floor(POSITION_X / bin_size) * bin_size,
           y_bin = floor(POSITION_Y / bin_size) * bin_size) %>%
    group_by(x_bin, y_bin, software) %>%
    summarise(mean_speed = mean(inst_speed, na.rm = TRUE),
              n = n(), .groups = "drop")
}

heat_data <- make_heatmap_data(all_vel, BIN_SIZE)

# Shared color scale
speed_range <- range(heat_data$mean_speed, na.rm = TRUE)

p_heat_tm <- ggplot(heat_data %>% filter(software == "TrackMate"),
                    aes(x = x_bin, y = y_bin, fill = mean_speed)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Speed\n(µm/min)", option = "inferno",
                       limits = speed_range) +
  coord_fixed() +
  labs(title = "TrackMate", x = "X (µm)", y = "Y (µm)") +
  theme_compare()

p_heat_im <- ggplot(heat_data %>% filter(software == "Imaris"),
                    aes(x = x_bin, y = y_bin, fill = mean_speed)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Speed\n(µm/min)", option = "inferno",
                       limits = speed_range) +
  coord_fixed() +
  labs(title = "Imaris", x = "X (µm)", y = "Y (µm)") +
  theme_compare()

heatmap_combined <- (p_heat_tm + p_heat_im) +
  plot_annotation(
    title = "Spatial Speed Heatmap Comparison (XY)",
    subtitle = sprintf("Bin size: %d µm | Shared color scale", BIN_SIZE),
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey50"))
  )

ggsave(file.path(OUTPUT_DIR, "09_spatial_speed_heatmaps.pdf"), heatmap_combined,
       width = 14, height = 7)
cat("  Saved: 09_spatial_speed_heatmaps.pdf\n")

# #############################################################################
# PART 12: INTERACTIVE 3D OVERLAY
# #############################################################################

cat("\n=== PART 12: INTERACTIVE 3D COMPARISON ===\n")

# Sample spots for interactive plot
MAX_INT <- 30000
set.seed(42)

tm_int <- tm_spots_clean %>%
  slice_sample(n = min(MAX_INT, nrow(tm_spots_clean))) %>%
  mutate(software = "TrackMate")
im_int <- im_spots_clean %>%
  slice_sample(n = min(MAX_INT, nrow(im_spots_clean))) %>%
  mutate(software = "Imaris")

fig_3d_overlay <- plot_ly() %>%
  add_trace(data = tm_int,
            x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
            type = "scatter3d", mode = "markers",
            marker = list(size = 1.5, opacity = 0.3, color = software_colors["TrackMate"]),
            name = "TrackMate") %>%
  add_trace(data = im_int,
            x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
            type = "scatter3d", mode = "markers",
            marker = list(size = 1.5, opacity = 0.3, color = software_colors["Imaris"]),
            name = "Imaris") %>%
  plotly::layout(
    title = list(text = "3D Spot Overlay: TrackMate (blue) vs Imaris (red)",
                 font = list(size = 16)),
    scene = list(
      xaxis = list(title = "X (µm)"),
      yaxis = list(title = "Y (µm)"),
      zaxis = list(title = "Z (µm)"),
      aspectmode = "data"
    )
  )

htmlwidgets::saveWidget(fig_3d_overlay,
                        file.path(OUTPUT_DIR, "10_interactive_3d_overlay.html"),
                        selfcontained = TRUE)
cat("  Saved: 10_interactive_3d_overlay.html\n")

# #############################################################################
# PART 13: SUMMARY TABLE & EXPORT
# #############################################################################

cat("\n=== PART 13: SUMMARY & EXPORT ===\n")

# Build comprehensive comparison table
comparison_table <- tibble(
  Metric = c(
    "Total spots (raw)",
    "Total tracks (raw)",
    "Tracks after QC",
    "Spots after QC",
    "Median spots/track",
    "Median duration (min)",
    "Median speed (µm/min)",
    "Mean speed (µm/min)",
    "Median displacement (µm)",
    "Median confinement ratio",
    "Median path length (µm)",
    "MSD α exponent",
    "% Directed",
    "% Random",
    "% Confined"
  ),
  TrackMate = c(
    nrow(tm_spots),
    nrow(tm_tracks),
    nrow(tm_tracks_clean),
    nrow(tm_spots_clean),
    median(tm_tracks_clean$N_SPOTS, na.rm = TRUE),
    median(tm_tracks_clean$DURATION_MIN, na.rm = TRUE),
    median(tm_tracks_clean$SPEED_MEAN_UM_MIN, na.rm = TRUE),
    mean(tm_tracks_clean$SPEED_MEAN_UM_MIN, na.rm = TRUE),
    median(tm_tracks_clean$TRACK_DISPLACEMENT, na.rm = TRUE),
    median(tm_tracks_clean$CONFINEMENT_RATIO, na.rm = TRUE),
    median(tm_tracks_clean$TOTAL_DISTANCE_TRAVELED, na.rm = TRUE),
    round(tm_msd$alpha, 3),
    move_summary %>% filter(software == "TrackMate", movement_type == "Directed") %>% pull(pct),
    move_summary %>% filter(software == "TrackMate", movement_type == "Random") %>% pull(pct),
    move_summary %>% filter(software == "TrackMate", movement_type == "Confined") %>% pull(pct)
  ),
  Imaris = c(
    nrow(im_spots),
    nrow(im_tracks),
    nrow(im_tracks_clean),
    nrow(im_spots_clean),
    median(im_tracks_clean$N_SPOTS, na.rm = TRUE),
    median(im_tracks_clean$DURATION_MIN, na.rm = TRUE),
    median(im_tracks_clean$SPEED_MEAN_UM_MIN, na.rm = TRUE),
    mean(im_tracks_clean$SPEED_MEAN_UM_MIN, na.rm = TRUE),
    median(im_tracks_clean$TRACK_DISPLACEMENT, na.rm = TRUE),
    median(im_tracks_clean$CONFINEMENT_RATIO, na.rm = TRUE),
    median(im_tracks_clean$TRACK_LENGTH, na.rm = TRUE),
    round(im_msd$alpha, 3),
    move_summary %>% filter(software == "Imaris", movement_type == "Directed") %>% pull(pct),
    move_summary %>% filter(software == "Imaris", movement_type == "Random") %>% pull(pct),
    move_summary %>% filter(software == "Imaris", movement_type == "Confined") %>% pull(pct)
  )
) %>%
  mutate(
    Ratio_TM_IM = round(TrackMate / Imaris, 3),
    Diff_pct    = round(100 * (TrackMate - Imaris) / pmax(TrackMate, Imaris), 1)
  )

cat("\n")
cat(strrep("=", 80), "\n")
cat("  TRACKMATE vs IMARIS — COMPARISON SUMMARY\n")
cat(strrep("=", 80), "\n\n")

print(comparison_table, n = Inf, width = Inf)

# Export
write_csv(comparison_table, file.path(OUTPUT_DIR, "comparison_summary.csv"))
write_csv(track_stats_table, file.path(OUTPUT_DIR, "track_stats_by_software.csv"))
write_csv(move_summary, file.path(OUTPUT_DIR, "movement_types_by_software.csv"))
write_csv(match_summary, file.path(OUTPUT_DIR, "spot_matching_summary.csv"))

cat("\n  Exported CSV files to", OUTPUT_DIR, "\n")

# #############################################################################
# FINAL OUTPUT LIST
# #############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("  OUTPUT FILES (in ", OUTPUT_DIR, "/)\n")
cat(strrep("=", 80), "\n\n")
cat("  PDFs:\n")
cat("    01_detection_comparison.pdf       — Spots per frame, XYZ distributions\n")
cat("    02_spot_matching.pdf              — Nearest-neighbor spatial matching\n")
cat("    03_track_distributions.pdf        — N spots, duration, speed, displacement, etc.\n")
cat("    04_velocity_comparison.pdf        — Instantaneous speed, turning angles\n")
cat("    05_msd_comparison.pdf             — Mean squared displacement overlay\n")
cat("    06_movement_classification.pdf    — Directed/Random/Confined proportions\n")
cat("    07_spatial_tracks.pdf             — Side-by-side track projections\n")
cat("    08_z_depth_comparison.pdf         — Speed by Z-layer\n")
cat("    09_spatial_speed_heatmaps.pdf     — XY speed heatmaps\n")
cat("\n  Interactive HTML:\n")
cat("    10_interactive_3d_overlay.html    — 3D spot overlay\n")
cat("\n  Data:\n")
cat("    comparison_summary.csv            — Full metric comparison table\n")
cat("    track_stats_by_software.csv       — Track statistics\n")
cat("    movement_types_by_software.csv    — Movement classification\n")
cat("    spot_matching_summary.csv         — Nearest-neighbor matching stats\n")

cat("\n")
cat(strrep("=", 80), "\n")
cat("  COMPARISON COMPLETE\n")
cat(strrep("=", 80), "\n")
