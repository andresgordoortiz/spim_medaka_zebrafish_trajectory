# =============================================================================
# TrackMate Filter & Validate
# =============================================================================
#
# PURPOSE:
#   This is Step 2 of the production pipeline:
#     Step 0: TrackMate → export ALL tracks (no filters)
#     Step 1: orientation_interactive.R → orient embryo → oriented CSVs
#   → Step 2: THIS SCRIPT → filter + validate → filtered CSVs
#     Step 3: trackmate_analysis.R → full biological analysis
#
# WHAT IT DOES:
#   1. Loads the oriented TrackMate CSVs (already reoriented)
#   2. Loads manual ground truth (40 cells + 5 YSL), applies drift correction
#   3. Loads Imaris data for 3-way comparison
#   4. Sweeps over QC filter parameters to find the combination that best
#      matches ground truth speed metrics
#   5. Produces diagnostic plots and a 3-way comparison
#   6. Exports filtered_tracks.csv and filtered_spots.csv for trackmate_analysis.R
#
# INPUTS:
#   - oriented_spots.csv  (from orientation_interactive.R)
#   - oriented_tracks.csv (from orientation_interactive.R)
#       OR 4D_spots.csv / 4D_tracks.csv if orientation was skipped
#   - high_res_mannualtracking/  (manual ground truth)
#   - all_tracks/  (Imaris tracking results)
#
# OUTPUTS (all in analysis_output/):
#   - 00_ground_truth_drift.pdf
#   - 00_ground_truth_summary.pdf
#   - 01_filter_impact.pdf
#   - 02_parameter_heatmaps.pdf
#   - 03_pareto_front.pdf
#   - 04_speed_comparison.pdf
#   - 05_filtering_effect.pdf
#   - filtered_tracks.csv          → input for trackmate_analysis.R
#   - filtered_spots.csv           → input for trackmate_analysis.R
#   - sweep_results.csv
#   - recommended_params.csv
#   - ground_truth_summary.csv
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
library(scales)

# =============================================================================
# PARAMETERS — edit these to match your experiment
# =============================================================================

FRAME_INTERVAL_SEC <- 30
FRAME_INTERVAL_MIN <- FRAME_INTERVAL_SEC / 60
SPEED_CONVERSION   <- 60 / FRAME_INTERVAL_SEC   # µm/frame → µm/min

MANUAL_DIR <- "high_res_mannualtracking"
IMARIS_DIR <- "all_tracks"
OUTPUT_DIR <- "analysis_output"

# Input files: use oriented CSVs if available, otherwise raw TrackMate
SPOTS_FILE  <- if (file.exists("oriented_spots.csv")) "oriented_spots.csv" else "4D_spots.csv"
TRACKS_FILE <- if (file.exists("oriented_tracks.csv")) "oriented_tracks.csv" else "4D_tracks.csv"

# --- Selection method ---
#   "pareto" = balanced trade-off (Pareto knee)
#   "score"  = prioritise ground truth similarity, then max tracks
SELECTION_METHOD <- "score"
SCORE_TOLERANCE  <- 0.05     # only for "score" mode: 0.05 = within 5% of GT speeds

# --- Sweep grid ---
SWEEP_GRID <- list(
  min_spots       = c(3, 5, 10, 15, 20, 30, 50, 75, 100),
  max_mean_speed  = c(3, 4, 5, 7, 10, Inf),
  max_max_speed   = c(5, 8, 10, 15, Inf),
  max_gaps        = c(0, 1, 2, 5, Inf),
  quality_pctile  = c(0, 0.1, 0.2, 0.3)
)

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# =============================================================================
# THEME
# =============================================================================

software_colors <- c("Manual (GT)"  = "#4DAF4A",
                     "TM optimized" = "#2166AC",
                     "TM default"   = "grey50",
                     "Imaris"       = "#B2182B")

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

cat("\n")
cat(strrep("=", 70), "\n")
cat("  TRACKMATE FILTER & VALIDATE\n")
cat(strrep("=", 70), "\n")
cat(sprintf("  Input:  %s / %s\n", SPOTS_FILE, TRACKS_FILE))
cat(sprintf("  Frame interval: %d sec | Speed conversion: ×%.0f\n",
            FRAME_INTERVAL_SEC, SPEED_CONVERSION))

# #############################################################################
# STEP 1: LOAD TRACKMATE DATA
# #############################################################################

cat("\n=== STEP 1: LOADING TRACKMATE DATA ===\n")

# Detect if files are oriented (have extra columns or standard TM format)
is_oriented <- grepl("oriented", SPOTS_FILE)

if (is_oriented) {
  # Oriented CSVs don't have the 3 metadata rows
  tm_spots <- read_csv(SPOTS_FILE, show_col_types = FALSE) %>%
    mutate(across(c(ID, TRACK_ID, POSITION_X, POSITION_Y, POSITION_Z, FRAME),
                  ~as.numeric(.))) %>%
    filter(!is.na(TRACK_ID))

  tm_tracks <- read_csv(TRACKS_FILE, show_col_types = FALSE) %>%
    mutate(across(-LABEL, ~suppressWarnings(as.numeric(.))))
} else {
  # Raw TrackMate CSVs: skip 3 metadata rows
  tm_spots <- read_csv(SPOTS_FILE, show_col_types = FALSE)[-c(1:3), ] %>%
    mutate(across(c(ID, TRACK_ID, POSITION_X, POSITION_Y, POSITION_Z, FRAME),
                  ~as.numeric(.))) %>%
    filter(!is.na(TRACK_ID))

  tm_tracks <- read_csv(TRACKS_FILE, show_col_types = FALSE)[-c(1:3), ] %>%
    mutate(across(-LABEL, ~suppressWarnings(as.numeric(.))))
}

# Convert units
tm_tracks <- tm_tracks %>%
  mutate(
    DURATION_MIN      = TRACK_DURATION * FRAME_INTERVAL_MIN,
    SPEED_MEAN_UM_MIN = TRACK_MEAN_SPEED * SPEED_CONVERSION,
    SPEED_MAX_UM_MIN  = TRACK_MAX_SPEED  * SPEED_CONVERSION,
    SPEED_STD_UM_MIN  = TRACK_STD_SPEED  * SPEED_CONVERSION,
    N_SPOTS           = NUMBER_SPOTS
  )

# Compute instantaneous speeds (for all spots)
tm_all_vel <- tm_spots %>%
  arrange(TRACK_ID, FRAME) %>%
  group_by(TRACK_ID) %>%
  mutate(
    dx = POSITION_X - lag(POSITION_X),
    dy = POSITION_Y - lag(POSITION_Y),
    dz = POSITION_Z - lag(POSITION_Z),
    dt_min = (FRAME - lag(FRAME)) * FRAME_INTERVAL_MIN,
    inst_speed = sqrt(dx^2 + dy^2 + dz^2) / dt_min
  ) %>%
  ungroup() %>%
  filter(!is.na(inst_speed) & dt_min > 0)

cat(sprintf("  Loaded %s tracks, %s spots\n",
            format(nrow(tm_tracks), big.mark = ","),
            format(nrow(tm_spots), big.mark = ",")))
cat(sprintf("  Unfiltered inst. speed median: %.3f µm/min\n",
            median(tm_all_vel$inst_speed, na.rm = TRUE)))

# #############################################################################
# STEP 2: LOAD & PROCESS MANUAL GROUND TRUTH
# #############################################################################

cat("\n=== STEP 2: MANUAL GROUND TRUTH (YSL DRIFT CORRECTION) ===\n")

# Helper: read one Imaris Position CSV
read_manual_position <- function(folder_name) {
  prefix <- sub("_Statistics$", "", folder_name)
  path   <- file.path(MANUAL_DIR, folder_name, paste0(prefix, "_Position.csv"))
  if (!file.exists(path)) { warning(sprintf("Missing: %s", path)); return(NULL) }

  read_csv(path, skip = 3, show_col_types = FALSE) %>%
    select(where(~ !all(is.na(.x)))) %>%
    rename(POSITION_X = `Position X`, POSITION_Y = `Position Y`,
           POSITION_Z = `Position Z`, FRAME = Time) %>%
    mutate(across(c(POSITION_X, POSITION_Y, POSITION_Z, FRAME), as.numeric),
           cell_name = prefix) %>%
    select(POSITION_X, POSITION_Y, POSITION_Z, FRAME, cell_name)
}

# Discover & load
all_folders <- list.dirs(MANUAL_DIR, recursive = FALSE, full.names = FALSE)
cell_folders <- grep("^Spots_cell_", all_folders, value = TRUE)
ysl_folders  <- grep("^Spots_nuclei_ysl_", all_folders, value = TRUE)

cell_data <- map_dfr(cell_folders, read_manual_position) %>% mutate(cell_type = "cell")
ysl_data  <- map_dfr(ysl_folders, read_manual_position)  %>% mutate(cell_type = "ysl")

cat(sprintf("  Loaded %d cells, %d YSL tracks\n",
            n_distinct(cell_data$cell_name), n_distinct(ysl_data$cell_name)))

# --- YSL drift correction ---
ysl_frame_disp <- ysl_data %>%
  arrange(cell_name, FRAME) %>%
  group_by(cell_name) %>%
  mutate(dx = POSITION_X - lag(POSITION_X),
         dy = POSITION_Y - lag(POSITION_Y),
         dz = POSITION_Z - lag(POSITION_Z)) %>%
  ungroup() %>%
  filter(!is.na(dx))

drift_per_frame <- ysl_frame_disp %>%
  group_by(FRAME) %>%
  summarise(drift_dx = mean(dx), drift_dy = mean(dy), drift_dz = mean(dz),
            drift_sd_x = sd(dx), drift_sd_y = sd(dy), drift_sd_z = sd(dz),
            n_ysl = n(), .groups = "drop") %>%
  mutate(drift_magnitude = sqrt(drift_dx^2 + drift_dy^2 + drift_dz^2),
         time_min = FRAME * FRAME_INTERVAL_MIN)

cat(sprintf("  Mean drift: %.3f µm/frame (%.3f µm/min)\n",
            mean(drift_per_frame$drift_magnitude),
            mean(drift_per_frame$drift_magnitude) / FRAME_INTERVAL_MIN))

# Apply drift correction to cells
correct_cell <- function(df, drift_df) {
  df %>%
    arrange(FRAME) %>%
    mutate(raw_dx = POSITION_X - lag(POSITION_X),
           raw_dy = POSITION_Y - lag(POSITION_Y),
           raw_dz = POSITION_Z - lag(POSITION_Z),
           dt_frames = FRAME - lag(FRAME)) %>%
    filter(!is.na(raw_dx) & dt_frames == 1) %>%
    left_join(drift_df %>% select(FRAME, drift_dx, drift_dy, drift_dz), by = "FRAME") %>%
    mutate(corr_dx = raw_dx - drift_dx,
           corr_dy = raw_dy - drift_dy,
           corr_dz = raw_dz - drift_dz,
           raw_speed  = sqrt(raw_dx^2 + raw_dy^2 + raw_dz^2) / FRAME_INTERVAL_MIN,
           corr_speed = sqrt(corr_dx^2 + corr_dy^2 + corr_dz^2) / FRAME_INTERVAL_MIN,
           time_min = FRAME * FRAME_INTERVAL_MIN)
}

cells_corrected <- cell_data %>%
  group_by(cell_name) %>%
  group_modify(~ correct_cell(.x, drift_per_frame)) %>%
  ungroup()

ysl_corrected <- ysl_data %>%
  group_by(cell_name) %>%
  group_modify(~ correct_cell(.x, drift_per_frame)) %>%
  ungroup()

# Per-cell summaries
cell_summaries <- cells_corrected %>%
  group_by(cell_name) %>%
  summarise(
    n_steps = n(), duration_min = n_steps * FRAME_INTERVAL_MIN,
    mean_corr_speed = mean(corr_speed, na.rm = TRUE),
    median_corr_speed = median(corr_speed, na.rm = TRUE),
    total_corr_dist = sum(sqrt(corr_dx^2 + corr_dy^2 + corr_dz^2), na.rm = TRUE),
    net_corr_disp = sqrt(sum(corr_dx, na.rm = TRUE)^2 +
                         sum(corr_dy, na.rm = TRUE)^2 +
                         sum(corr_dz, na.rm = TRUE)^2),
    confinement_corr = net_corr_disp / (total_corr_dist + 1e-10),
    .groups = "drop")

# Ground truth targets
GT <- list(
  inst_speed_median      = median(cells_corrected$corr_speed, na.rm = TRUE),
  inst_speed_mean        = mean(cells_corrected$corr_speed, na.rm = TRUE),
  track_mean_speed_median = median(cell_summaries$mean_corr_speed, na.rm = TRUE)
)

cat(sprintf("\n  Ground truth targets (drift-corrected, %d cells):\n", nrow(cell_summaries)))
cat(sprintf("    Inst. speed median:      %.3f µm/min\n", GT$inst_speed_median))
cat(sprintf("    Inst. speed mean:        %.3f µm/min\n", GT$inst_speed_mean))
cat(sprintf("    Per-track speed median:  %.3f µm/min\n", GT$track_mean_speed_median))
cat(sprintf("    YSL corrected speed:     %.3f µm/min (should be ≈0)\n",
            median(ysl_corrected$corr_speed, na.rm = TRUE)))

# --- Ground truth diagnostic plots ---

# Drift correction
manual_all <- bind_rows(cell_data, ysl_data)

p_raw_xy <- ggplot(manual_all, aes(x = POSITION_X, y = POSITION_Y,
                                    group = cell_name, color = cell_type)) +
  geom_path(linewidth = 0.5, alpha = 0.7) +
  scale_color_manual(values = c("cell" = "#2166AC", "ysl" = "#E66101"),
                     labels = c("Cell", "YSL (reference)")) +
  coord_fixed() +
  labs(title = "Raw XY Tracks", subtitle = "YSL apparent motion = embryo drift",
       x = "X (µm)", y = "Y (µm)", color = NULL) +
  theme_pipe()

cum_drift <- drift_per_frame %>%
  arrange(FRAME) %>%
  mutate(cum_x = cumsum(drift_dx), cum_y = cumsum(drift_dy))

p_drift_path <- ggplot(cum_drift, aes(x = cum_x, y = cum_y, color = time_min)) +
  geom_path(linewidth = 1.2) + geom_point(size = 1) +
  scale_color_viridis_c(name = "Time\n(min)") + coord_fixed() +
  labs(title = "Cumulative Embryo Drift (XY)", x = "ΔX (µm)", y = "ΔY (µm)") +
  theme_pipe()

drift_long <- drift_per_frame %>%
  pivot_longer(cols = c(drift_dx, drift_dy, drift_dz),
               names_to = "axis", values_to = "drift") %>%
  mutate(axis = sub("drift_d", "", axis) %>% toupper())

p_drift_axes <- ggplot(drift_long, aes(x = time_min, y = drift, color = axis)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("X" = "#E41A1C", "Y" = "#377EB8", "Z" = "#4DAF4A")) +
  labs(title = "Drift Per Axis", x = "Time (min)", y = "Drift (µm)", color = "Axis") +
  theme_pipe()

speed_compare <- bind_rows(
  cells_corrected %>% transmute(speed = raw_speed, type = "Raw (with drift)"),
  cells_corrected %>% transmute(speed = corr_speed, type = "Corrected")
)

p_speed_corr <- ggplot(speed_compare, aes(x = speed, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Raw (with drift)" = "grey60", "Corrected" = "#4DAF4A")) +
  coord_cartesian(xlim = c(0, quantile(speed_compare$speed, 0.99))) +
  labs(title = "Raw vs Corrected Speed", x = "Speed (µm/min)", y = "Density", fill = NULL) +
  theme_pipe()

p_gt_drift <- (p_raw_xy + p_drift_path) / (p_drift_axes + p_speed_corr) +
  plot_annotation(title = "Ground Truth: YSL Drift Correction",
                  theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))

ggsave(file.path(OUTPUT_DIR, "00_ground_truth_drift.pdf"), p_gt_drift, width = 14, height = 10)
cat("  Saved: 00_ground_truth_drift.pdf\n")

# Per-cell summary plot
p_cell_speeds <- ggplot(cell_summaries %>% arrange(mean_corr_speed) %>%
                          mutate(cell_name = factor(cell_name, levels = cell_name)),
                        aes(x = cell_name, y = mean_corr_speed)) +
  geom_col(fill = "#4DAF4A", alpha = 0.8) + coord_flip() +
  labs(title = "Per-Cell Mean Speed (corrected)", x = NULL, y = "µm/min") +
  theme_pipe() + theme(axis.text.y = element_text(size = 6))

p_confine <- ggplot(cell_summaries, aes(x = confinement_corr)) +
  geom_histogram(bins = 15, fill = "#4DAF4A", alpha = 0.7) +
  geom_vline(xintercept = c(0.3, 0.7), linetype = "dashed", color = "grey50") +
  labs(title = "Confinement Ratio (corrected)", x = "Confinement", y = "Count") +
  theme_pipe()

p_gt_summary <- (p_cell_speeds + p_confine) +
  plot_annotation(title = "Ground Truth: Per-Cell Summary",
                  theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))

ggsave(file.path(OUTPUT_DIR, "00_ground_truth_summary.pdf"), p_gt_summary, width = 14, height = 6)
cat("  Saved: 00_ground_truth_summary.pdf\n")

# #############################################################################
# STEP 3: LOAD IMARIS (for comparison)
# #############################################################################

cat("\n=== STEP 3: LOADING IMARIS DATA ===\n")

imaris_speeds <- NULL
imaris_track_speeds <- NULL
imaris_n_tracks <- 0

FILE_PREFIX <- "ome-tiff.companion-bin-8bit-crop"

read_imaris_csv <- function(stat_name) {
  path <- file.path(IMARIS_DIR, paste0(FILE_PREFIX, "_", stat_name, ".csv"))
  if (!file.exists(path)) return(NULL)
  read_csv(path, skip = 3, show_col_types = FALSE) %>%
    select(where(~ !all(is.na(.x))))
}

im_pos <- read_imaris_csv("Position")

if (!is.null(im_pos)) {
  im_spots <- im_pos %>%
    rename(POSITION_X = `Position X`, POSITION_Y = `Position Y`,
           POSITION_Z = `Position Z`, FRAME = Time,
           TRACK_ID = TrackID) %>%
    mutate(across(c(TRACK_ID, POSITION_X, POSITION_Y, POSITION_Z, FRAME), as.numeric))

  im_spots_per_track <- im_spots %>% count(TRACK_ID, name = "N_SPOTS")

  im_track_spd <- read_imaris_csv("Track_Speed_Mean")
  if (!is.null(im_track_spd)) {
    im_tracks <- im_track_spd %>%
      rename(TRACK_SPEED_MEAN_RAW = `Track Speed Mean`, TRACK_ID = ID) %>%
      mutate(across(everything(), as.numeric)) %>%
      left_join(im_spots_per_track, by = "TRACK_ID") %>%
      mutate(SPEED_MEAN_UM_MIN = TRACK_SPEED_MEAN_RAW / FRAME_INTERVAL_SEC * 60) %>%
      filter(N_SPOTS >= 5 & SPEED_MEAN_UM_MIN <= 10)

    imaris_track_speeds <- im_tracks$SPEED_MEAN_UM_MIN
    imaris_n_tracks <- nrow(im_tracks)

    im_vel <- im_spots %>%
      filter(TRACK_ID %in% im_tracks$TRACK_ID) %>%
      arrange(TRACK_ID, FRAME) %>%
      group_by(TRACK_ID) %>%
      mutate(dx = POSITION_X - lag(POSITION_X),
             dy = POSITION_Y - lag(POSITION_Y),
             dz = POSITION_Z - lag(POSITION_Z),
             dt_min = (FRAME - lag(FRAME)) * FRAME_INTERVAL_MIN,
             inst_speed = sqrt(dx^2 + dy^2 + dz^2) / dt_min) %>%
      ungroup() %>%
      filter(!is.na(inst_speed) & dt_min > 0)

    imaris_speeds <- im_vel$inst_speed
  }

  cat(sprintf("  Imaris: %s tracks (after QC ≥5 spots, ≤10 µm/min)\n",
              format(imaris_n_tracks, big.mark = ",")))
} else {
  cat("  Imaris data not found — skipping comparison\n")
}

# #############################################################################
# STEP 4: PARAMETER SWEEP
# #############################################################################

cat("\n=== STEP 4: PARAMETER SWEEP ===\n")

param_grid <- expand.grid(SWEEP_GRID, stringsAsFactors = FALSE)
cat(sprintf("  Sweeping %s filter combinations...\n",
            format(nrow(param_grid), big.mark = ",")))

# Pre-split spots by track for fast lookup
vel_by_track <- split(tm_all_vel$inst_speed, tm_all_vel$TRACK_ID)

# Evaluation function
eval_combo <- function(min_spots, max_mean_speed, max_max_speed, max_gaps, quality_pctile) {
  q_thresh <- quantile(tm_tracks$TRACK_MEAN_QUALITY, quality_pctile, na.rm = TRUE)

  passing <- tm_tracks %>%
    filter(NUMBER_SPOTS >= min_spots,
           SPEED_MEAN_UM_MIN <= max_mean_speed,
           SPEED_MAX_UM_MIN  <= max_max_speed,
           NUMBER_GAPS <= max_gaps,
           TRACK_MEAN_QUALITY >= q_thresh)

  n <- nrow(passing)
  if (n < 5) return(list(n_tracks = n, inst_speed_median = NA,
                          inst_speed_mean = NA, track_speed_median = NA))

  passing_ids <- as.character(passing$TRACK_ID)
  all_speeds <- unlist(vel_by_track[passing_ids], use.names = FALSE)

  list(n_tracks = n,
       inst_speed_median  = median(all_speeds, na.rm = TRUE),
       inst_speed_mean    = mean(all_speeds, na.rm = TRUE),
       track_speed_median = median(passing$SPEED_MEAN_UM_MIN, na.rm = TRUE))
}

# Run sweep
sweep_results <- tibble(
  min_spots      = param_grid$min_spots,
  max_mean_speed = param_grid$max_mean_speed,
  max_max_speed  = param_grid$max_max_speed,
  max_gaps       = param_grid$max_gaps,
  quality_pctile = param_grid$quality_pctile,
  n_tracks = NA_real_, inst_speed_median = NA_real_,
  inst_speed_mean = NA_real_, track_speed_median = NA_real_
)

pb <- txtProgressBar(min = 0, max = nrow(param_grid), style = 3)
for (i in seq_len(nrow(param_grid))) {
  r <- eval_combo(param_grid$min_spots[i], param_grid$max_mean_speed[i],
                  param_grid$max_max_speed[i], param_grid$max_gaps[i],
                  param_grid$quality_pctile[i])
  sweep_results$n_tracks[i]          <- r$n_tracks
  sweep_results$inst_speed_median[i] <- r$inst_speed_median
  sweep_results$inst_speed_mean[i]   <- r$inst_speed_mean
  sweep_results$track_speed_median[i] <- r$track_speed_median
  setTxtProgressBar(pb, i)
}
close(pb)

cat(sprintf("  Done. Valid combos: %d\n",
            sum(!is.na(sweep_results$inst_speed_median))))

# #############################################################################
# STEP 5: SCORING & SELECTION
# #############################################################################

cat("\n=== STEP 5: SCORING & SELECTION ===\n")
cat(sprintf("  Method: %s\n", SELECTION_METHOD))

N_TOTAL_TRACKS <- nrow(tm_tracks)

sweep_scored <- sweep_results %>%
  filter(!is.na(inst_speed_median)) %>%
  mutate(
    err_inst_median = abs(inst_speed_median - GT$inst_speed_median) / GT$inst_speed_median,
    err_inst_mean   = abs(inst_speed_mean   - GT$inst_speed_mean)   / GT$inst_speed_mean,
    err_track_speed = abs(track_speed_median - GT$track_mean_speed_median) / GT$track_mean_speed_median,
    score = 0.4 * err_inst_median + 0.3 * err_inst_mean + 0.3 * err_track_speed,
    frac_tracks = n_tracks / N_TOTAL_TRACKS
  ) %>%
  arrange(score)

# Pareto front (always computed for visualisation)
identify_pareto <- function(df) {
  df <- df %>% arrange(score)
  is_pareto <- rep(TRUE, nrow(df))
  max_tracks <- 0
  for (i in seq_len(nrow(df))) {
    if (df$n_tracks[i] >= max_tracks) { max_tracks <- df$n_tracks[i] }
    else { is_pareto[i] <- FALSE }
  }
  df$is_pareto <- is_pareto
  df
}

sweep_scored <- identify_pareto(sweep_scored)
pareto_front <- sweep_scored %>% filter(is_pareto) %>% arrange(score)

cat(sprintf("  Valid combos: %d | Pareto-optimal: %d\n",
            nrow(sweep_scored), nrow(pareto_front)))

# Select best
if (SELECTION_METHOD == "pareto") {
  pareto_norm <- pareto_front %>%
    mutate(score_norm  = (score - min(score)) / (max(score) - min(score) + 1e-9),
           tracks_norm = (n_tracks - min(n_tracks)) / (max(n_tracks) - min(n_tracks) + 1e-9),
           dist_to_ideal = sqrt(score_norm^2 + (1 - tracks_norm)^2))
  best_idx <- which.min(pareto_norm$dist_to_ideal)
  best <- pareto_norm[best_idx, ]
  score_cutoff <- NA
  cat(sprintf("  Pareto knee: score=%.4f, %s tracks\n",
              best$score, format(best$n_tracks, big.mark = ",")))
} else {
  score_cutoff <- SCORE_TOLERANCE
  near_best <- sweep_scored %>% filter(score <= score_cutoff) %>% arrange(desc(n_tracks))
  if (nrow(near_best) == 0) {
    cat(sprintf("  WARNING: No combos within %.0f%% of GT. Using best available.\n",
                SCORE_TOLERANCE * 100))
    best <- sweep_scored[1, ]
    score_cutoff <- best$score
  } else {
    best <- near_best[1, ]
  }
  best_idx <- NA
  cat(sprintf("  GT tolerance: ±%.0f%% (score ≤ %.4f)\n", SCORE_TOLERANCE * 100, score_cutoff))
  cat(sprintf("  Combos in range: %d\n", nrow(near_best)))
  cat(sprintf("  Selected: score=%.4f, %s tracks\n",
              best$score, format(best$n_tracks, big.mark = ",")))
}

cat(sprintf("\n  RECOMMENDED PARAMETERS:\n"))
cat(sprintf("  ┌─────────────────────────────────────────────────────────────┐\n"))
cat(sprintf("  │  min_spots       = %-5d  (min detections per track)       │\n", best$min_spots))
cat(sprintf("  │  max_mean_speed  = %-5.1f  µm/min                         │\n", best$max_mean_speed))
cat(sprintf("  │  max_max_speed   = %-5.1f  µm/min                         │\n", best$max_max_speed))
cat(sprintf("  │  max_gaps        = %-5.0f  (max frame gaps)               │\n", best$max_gaps))
cat(sprintf("  │  quality_pctile  = %-5.0f%% (bottom quality removed)       │\n", best$quality_pctile * 100))
cat(sprintf("  └─────────────────────────────────────────────────────────────┘\n"))
cat(sprintf("\n  Closeness to ground truth:\n"))
cat(sprintf("    %-22s  %8s  %8s  %6s\n", "", "Filtered", "GT", "Error"))
cat(sprintf("    Inst speed median:  %8.3f  %8.3f  %5.1f%%\n",
            best$inst_speed_median, GT$inst_speed_median, best$err_inst_median * 100))
cat(sprintf("    Inst speed mean:    %8.3f  %8.3f  %5.1f%%\n",
            best$inst_speed_mean, GT$inst_speed_mean, best$err_inst_mean * 100))
cat(sprintf("    Track speed median: %8.3f  %8.3f  %5.1f%%\n",
            best$track_speed_median, GT$track_mean_speed_median, best$err_track_speed * 100))
cat(sprintf("  Tracks retained: %s / %s (%.1f%%)\n",
            format(best$n_tracks, big.mark = ","),
            format(N_TOTAL_TRACKS, big.mark = ","),
            best$frac_tracks * 100))

# #############################################################################
# STEP 6: APPLY BEST FILTERS
# #############################################################################

cat("\n=== STEP 6: APPLYING FILTERS ===\n")

q_thresh <- quantile(tm_tracks$TRACK_MEAN_QUALITY, best$quality_pctile, na.rm = TRUE)

tm_filtered <- tm_tracks %>%
  filter(NUMBER_SPOTS >= best$min_spots,
         SPEED_MEAN_UM_MIN <= best$max_mean_speed,
         SPEED_MAX_UM_MIN  <= best$max_max_speed,
         NUMBER_GAPS <= best$max_gaps,
         TRACK_MEAN_QUALITY >= q_thresh)

tm_filtered_spots <- tm_spots %>%
  filter(TRACK_ID %in% tm_filtered$TRACK_ID)

tm_filtered_vel <- tm_all_vel %>%
  filter(TRACK_ID %in% tm_filtered$TRACK_ID)

# Default comparison set (simple ≥5 spots, ≤10 µm/min)
tm_default <- tm_tracks %>% filter(NUMBER_SPOTS >= 5 & SPEED_MEAN_UM_MIN <= 10)
tm_default_vel <- tm_all_vel %>% filter(TRACK_ID %in% tm_default$TRACK_ID)

cat(sprintf("  Optimized: %s tracks, %s spots\n",
            format(nrow(tm_filtered), big.mark = ","),
            format(nrow(tm_filtered_spots), big.mark = ",")))
cat(sprintf("  Default:   %s tracks (for comparison)\n",
            format(nrow(tm_default), big.mark = ",")))

# #############################################################################
# STEP 7: DIAGNOSTIC PLOTS
# #############################################################################

cat("\n=== STEP 7: PLOTS ===\n")

# ─── PLOT 1: Filter impact (marginal effects — line plots) ───────────────────

# Compute medians for each parameter value, averaging over others
marginal_data <- bind_rows(
  sweep_scored %>%
    group_by(value = min_spots) %>%
    summarise(inst_median = median(inst_speed_median), track_median = median(track_speed_median),
              score_med = median(score), n_tracks_med = median(n_tracks), .groups = "drop") %>%
    mutate(param = "min_spots"),
  sweep_scored %>% filter(is.finite(max_mean_speed)) %>%
    group_by(value = max_mean_speed) %>%
    summarise(inst_median = median(inst_speed_median), track_median = median(track_speed_median),
              score_med = median(score), n_tracks_med = median(n_tracks), .groups = "drop") %>%
    mutate(param = "max_mean_speed"),
  sweep_scored %>% filter(is.finite(max_max_speed)) %>%
    group_by(value = max_max_speed) %>%
    summarise(inst_median = median(inst_speed_median), track_median = median(track_speed_median),
              score_med = median(score), n_tracks_med = median(n_tracks), .groups = "drop") %>%
    mutate(param = "max_max_speed"),
  sweep_scored %>% filter(is.finite(max_gaps)) %>%
    group_by(value = max_gaps) %>%
    summarise(inst_median = median(inst_speed_median), track_median = median(track_speed_median),
              score_med = median(score), n_tracks_med = median(n_tracks), .groups = "drop") %>%
    mutate(param = "max_gaps"),
  sweep_scored %>%
    group_by(value = quality_pctile * 100) %>%
    summarise(inst_median = median(inst_speed_median), track_median = median(track_speed_median),
              score_med = median(score), n_tracks_med = median(n_tracks), .groups = "drop") %>%
    mutate(param = "quality_%_removed")
) %>%
  mutate(param = factor(param, levels = c("min_spots", "max_mean_speed", "max_max_speed",
                                           "max_gaps", "quality_%_removed")))

# Row 1: Speed metrics vs GT
marginal_speed <- marginal_data %>%
  pivot_longer(cols = c(inst_median, track_median),
               names_to = "metric", values_to = "speed") %>%
  mutate(metric = recode(metric,
                          "inst_median"  = "Inst. speed median",
                          "track_median" = "Per-track speed median"))

p1_speed <- ggplot(marginal_speed, aes(x = value, y = speed, color = metric)) +
  geom_line(linewidth = 1) + geom_point(size = 1.5) +
  geom_hline(yintercept = GT$inst_speed_median, linetype = "dashed", color = "#4DAF4A", linewidth = 0.5) +
  geom_hline(yintercept = GT$track_mean_speed_median, linetype = "dotted", color = "#4DAF4A", linewidth = 0.5) +
  scale_color_manual(values = c("Inst. speed median" = "#2166AC", "Per-track speed median" = "#B2182B")) +
  facet_wrap(~param, scales = "free_x", nrow = 1) +
  labs(title = "How each filter parameter affects speed metrics",
       subtitle = "Green dashed/dotted = ground truth targets",
       x = "Parameter value", y = "Speed (µm/min)", color = NULL) +
  theme_pipe(base_size = 9) + theme(strip.text = element_text(size = 8))

# Row 2: Score (lower = better match to GT)
p1_score <- ggplot(marginal_data, aes(x = value, y = score_med)) +
  geom_line(linewidth = 0.8, color = "#762A83") + geom_point(size = 1.5, color = "#762A83") +
  facet_wrap(~param, scales = "free_x", nrow = 1) +
  labs(title = "Score (lower = closer to GT)",
       x = "Parameter value", y = "Score") +
  theme_pipe(base_size = 9) + theme(strip.text = element_text(size = 8))

# Row 3: Tracks retained
p1_n <- ggplot(marginal_data, aes(x = value, y = n_tracks_med)) +
  geom_line(linewidth = 0.8, color = "grey40") + geom_point(size = 1.5) +
  scale_y_log10() +
  facet_wrap(~param, scales = "free_x", nrow = 1) +
  labs(title = "Tracks remaining after each filter",
       x = "Parameter value", y = "N tracks (log scale)") +
  theme_pipe(base_size = 9) + theme(strip.text = element_text(size = 8))

p1 <- p1_speed / p1_score / p1_n +
  plot_annotation(
    title = "Filter Impact Analysis",
    subtitle = "Which QC parameters affect speed, score, and data retention?",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50"))
  )

ggsave(file.path(OUTPUT_DIR, "01_filter_impact.pdf"), p1, width = 16, height = 10)
cat("  Saved: 01_filter_impact.pdf\n")

# ─── PLOT 2: 2D Heatmaps (separate color scales) ─────────────────────────────

heatmap_data <- sweep_scored %>%
  filter(is.finite(max_mean_speed)) %>%
  group_by(min_spots, max_mean_speed) %>%
  summarise(mean_score = mean(score), mean_inst = mean(inst_speed_median),
            mean_tracks = mean(n_tracks), .groups = "drop")

p2a <- ggplot(heatmap_data, aes(x = factor(min_spots), y = factor(max_mean_speed),
                                 fill = mean_score)) +
  geom_tile() + geom_text(aes(label = sprintf("%.2f", mean_score)), size = 2.3, color = "white") +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Score") +
  labs(title = "Score (lower = closer to GT)", x = "Min spots", y = "Max mean speed") +
  theme_pipe(base_size = 9)

p2b <- ggplot(heatmap_data, aes(x = factor(min_spots), y = factor(max_mean_speed),
                                 fill = mean_inst)) +
  geom_tile() + geom_text(aes(label = sprintf("%.2f", mean_inst)), size = 2.3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = GT$inst_speed_median, name = "µm/min") +
  labs(title = sprintf("Inst. Speed Median (target: %.2f)", GT$inst_speed_median),
       x = "Min spots", y = "Max mean speed") +
  theme_pipe(base_size = 9)

p2c <- ggplot(heatmap_data, aes(x = factor(min_spots), y = factor(max_mean_speed),
                                 fill = mean_tracks)) +
  geom_tile() + geom_text(aes(label = format(round(mean_tracks), big.mark = ",")), size = 2.1) +
  scale_fill_viridis_c(name = "N tracks", labels = comma) +
  labs(title = "Tracks Retained", x = "Min spots", y = "Max mean speed") +
  theme_pipe(base_size = 9)

p2 <- p2a + p2b + p2c +
  plot_annotation(title = "Parameter Space: MIN_SPOTS × MAX_MEAN_SPEED",
                  subtitle = "Each heatmap has its own color scale",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                                plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50")))

ggsave(file.path(OUTPUT_DIR, "02_parameter_heatmaps.pdf"), p2, width = 18, height = 6)
cat("  Saved: 02_parameter_heatmaps.pdf\n")

# ─── PLOT 3: Pareto front ────────────────────────────────────────────────────

p3 <- ggplot(sweep_scored, aes(x = n_tracks, y = score)) +
  geom_point(alpha = 0.08, size = 0.8, color = "grey50") +
  geom_point(data = pareto_front, color = "#E41A1C", size = 2) +
  geom_line(data = pareto_front, color = "#E41A1C", linewidth = 0.8) +
  geom_point(data = best, aes(x = n_tracks, y = score),
             color = "#2166AC", size = 5, shape = 18) +
  {if (!is.na(score_cutoff))
    geom_hline(yintercept = score_cutoff, linetype = "dotted",
               color = "#2166AC", linewidth = 0.5)} +
  {if (!is.na(score_cutoff))
    annotate("text", x = min(sweep_scored$n_tracks), y = score_cutoff + 0.005,
             label = sprintf("within %.0f%% of GT speeds", SCORE_TOLERANCE * 100),
             hjust = 0, size = 2.5, color = "#2166AC")} +
  annotate("label", x = best$n_tracks, y = best$score + 0.02,
           label = sprintf("SELECTED\n%s tracks, score=%.3f",
                           format(best$n_tracks, big.mark = ","), best$score),
           size = 3, color = "#2166AC", fill = "white", label.size = 0.3) +
  scale_x_log10(labels = comma) +
  labs(title = "Pareto Front: GT Similarity vs Data Retention",
       subtitle = paste0("Red = Pareto front | Blue diamond = selected (",
                         if (SELECTION_METHOD == "pareto") "knee" else "within GT tolerance", ")"),
       x = "Tracks retained (log)", y = "Score (lower = closer to GT)") +
  theme_pipe()

ggsave(file.path(OUTPUT_DIR, "03_pareto_front.pdf"), p3, width = 10, height = 7)
cat("  Saved: 03_pareto_front.pdf\n")

# ─── PLOT 4: Speed comparison (faceted + overlay + box) ───────────────────────

speed_all <- bind_rows(
  tibble(speed = tm_default_vel$inst_speed, dataset = "TM default"),
  tibble(speed = tm_filtered_vel$inst_speed, dataset = "TM optimized"),
  if (!is.null(imaris_speeds)) tibble(speed = imaris_speeds, dataset = "Imaris") else NULL,
  tibble(speed = cells_corrected$corr_speed, dataset = "Manual (GT)")
) %>%
  mutate(dataset = factor(dataset, levels = c("TM default", "TM optimized", "Imaris", "Manual (GT)")))

speed_stats <- speed_all %>%
  group_by(dataset) %>%
  summarise(median = median(speed, na.rm = TRUE),
            mean = mean(speed, na.rm = TRUE),
            n = n(), .groups = "drop")

# Faceted histograms
p4a <- ggplot(speed_all, aes(x = speed, fill = dataset)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, alpha = 0.7) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  geom_vline(data = speed_stats, aes(xintercept = median),
             linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_vline(xintercept = GT$inst_speed_median,
             linetype = "dotted", color = "#4DAF4A", linewidth = 0.6) +
  geom_text(data = speed_stats,
            aes(x = Inf, y = Inf,
                label = sprintf("med=%.2f\nmean=%.2f\nn=%s", median, mean,
                                format(n, big.mark = ","))),
            hjust = 1.1, vjust = 1.3, size = 2.8, inherit.aes = FALSE) +
  facet_wrap(~dataset, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = software_colors) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(title = "Instantaneous Speed — Per Pipeline",
       subtitle = sprintf("Green dotted = GT target (%.2f µm/min)", GT$inst_speed_median),
       x = "Speed (µm/min)", y = "Density") +
  theme_pipe(base_size = 9) + theme(legend.position = "none")

# Overlay density
p4b <- ggplot(speed_all, aes(x = speed, color = dataset)) +
  geom_density(linewidth = 0.9) +
  geom_vline(xintercept = GT$inst_speed_median, linetype = "dashed",
             color = "#4DAF4A", linewidth = 0.5) +
  scale_color_manual(values = software_colors) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(title = "Overlay", x = "Speed (µm/min)", y = "Density", color = NULL) +
  theme_pipe()

# Box plot
p4c <- ggplot(speed_all, aes(x = dataset, y = speed, fill = dataset)) +
  geom_boxplot(alpha = 0.6, outlier.alpha = 0.02, outlier.size = 0.3) +
  geom_hline(yintercept = GT$inst_speed_median, linetype = "dashed",
             color = "#4DAF4A", linewidth = 0.6) +
  scale_fill_manual(values = software_colors) +
  coord_cartesian(ylim = c(0, 10)) +
  labs(title = "Box Plot", x = NULL, y = "Speed (µm/min)") +
  theme_pipe() + theme(legend.position = "none",
                       axis.text.x = element_text(angle = 15, hjust = 1))

p4 <- p4a / (p4b + p4c) +
  plot_annotation(title = "Speed Comparison Across Pipelines",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

ggsave(file.path(OUTPUT_DIR, "04_speed_comparison.pdf"), p4, width = 16, height = 10)
cat("  Saved: 04_speed_comparison.pdf\n")

# ─── PLOT 5: Filtering effect (duration + spots before/after) ────────────────

dur_compare <- bind_rows(
  tibble(duration = tm_tracks$DURATION_MIN, group = "All (unfiltered)"),
  tibble(duration = tm_default$DURATION_MIN, group = "TM default"),
  tibble(duration = tm_filtered$DURATION_MIN, group = "TM optimized")
) %>% mutate(group = factor(group, levels = c("All (unfiltered)", "TM default", "TM optimized")))

p5a <- ggplot(dur_compare, aes(x = duration, fill = group)) +
  geom_histogram(bins = 60, alpha = 0.7) +
  facet_wrap(~group, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("All (unfiltered)" = "grey80",
                                "TM default" = "grey50",
                                "TM optimized" = "#2166AC")) +
  labs(title = "Track Duration — Before vs After", x = "Duration (min)", y = "Count") +
  theme_pipe(base_size = 9) + theme(legend.position = "none")

spots_compare <- bind_rows(
  tibble(n_spots = tm_tracks$NUMBER_SPOTS, group = "All"),
  tibble(n_spots = tm_default$NUMBER_SPOTS, group = "Default"),
  tibble(n_spots = tm_filtered$NUMBER_SPOTS, group = "Optimized")
) %>% mutate(group = factor(group, levels = c("All", "Default", "Optimized")))

p5b <- ggplot(spots_compare, aes(x = n_spots, fill = group)) +
  geom_histogram(bins = 60, alpha = 0.7) +
  facet_wrap(~group, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("All" = "grey80", "Default" = "grey50",
                                "Optimized" = "#2166AC")) +
  labs(title = "Spots per Track", x = "N spots", y = "Count") +
  theme_pipe(base_size = 9) + theme(legend.position = "none")

p5 <- p5a / p5b +
  plot_annotation(title = "What the Filters Remove",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

ggsave(file.path(OUTPUT_DIR, "05_filtering_effect.pdf"), p5, width = 16, height = 8)
cat("  Saved: 05_filtering_effect.pdf\n")

# #############################################################################
# STEP 8: CONSOLE SUMMARY
# #############################################################################

cat("\n=== STEP 8: SUMMARY ===\n\n")

cat("  Speed comparison (inst. speed, µm/min):\n")
cat(sprintf("  %-18s  %8s  %8s  %10s\n", "Pipeline", "Median", "Mean", "N values"))
cat("  ", strrep("-", 50), "\n")
for (i in 1:nrow(speed_stats)) {
  r <- speed_stats[i, ]
  cat(sprintf("  %-18s  %8.3f  %8.3f  %10s\n",
              as.character(r$dataset), r$median, r$mean, format(r$n, big.mark = ",")))
}
cat(sprintf("  %-18s  %8.3f  %8.3f  %10s\n",
            "GROUND TRUTH", GT$inst_speed_median, GT$inst_speed_mean, "—"))

# #############################################################################
# STEP 9: EXPORT
# #############################################################################

cat("\n=== STEP 9: EXPORT ===\n")

# Filtered tracks
tracks_out <- tm_filtered %>%
  select(any_of(c("TRACK_ID", "LABEL", "NUMBER_SPOTS", "NUMBER_GAPS", "NUMBER_SPLITS",
                   "NUMBER_MERGES", "LONGEST_GAP", "TRACK_DURATION", "TRACK_START",
                   "TRACK_STOP", "TRACK_DISPLACEMENT", "TRACK_X_LOCATION",
                   "TRACK_Y_LOCATION", "TRACK_Z_LOCATION", "TRACK_MEAN_SPEED",
                   "TRACK_MAX_SPEED", "TRACK_MIN_SPEED", "TRACK_MEDIAN_SPEED",
                   "TRACK_STD_SPEED", "TRACK_MEAN_QUALITY", "TOTAL_DISTANCE_TRAVELED",
                   "MAX_DISTANCE_TRAVELED", "CONFINEMENT_RATIO",
                   "MEAN_STRAIGHT_LINE_SPEED", "LINEARITY_OF_FORWARD_PROGRESSION",
                   "MEAN_DIRECTIONAL_CHANGE_RATE",
                   "DURATION_MIN", "SPEED_MEAN_UM_MIN", "SPEED_MAX_UM_MIN", "SPEED_STD_UM_MIN")))

write_csv(tracks_out, file.path(OUTPUT_DIR, "filtered_tracks.csv"))

# Filtered spots
spots_out <- tm_filtered_spots %>%
  select(any_of(c("ID", "LABEL", "TRACK_ID", "QUALITY",
                   "POSITION_X", "POSITION_Y", "POSITION_Z",
                   "FRAME", "RADIUS", "VISIBILITY",
                   "MEAN_INTENSITY_CH1", "MEDIAN_INTENSITY_CH1",
                   "MIN_INTENSITY_CH1", "MAX_INTENSITY_CH1",
                   "TOTAL_INTENSITY_CH1", "STD_INTENSITY_CH1",
                   "CONTRAST_CH1", "SNR_CH1")))

write_csv(spots_out, file.path(OUTPUT_DIR, "filtered_spots.csv"))

# Sweep results + params
write_csv(sweep_scored, file.path(OUTPUT_DIR, "sweep_results.csv"))
write_csv(
  tibble(parameter   = c("min_spots", "max_mean_speed", "max_max_speed",
                          "max_gaps", "quality_pctile_removed"),
         value       = c(best$min_spots, best$max_mean_speed, best$max_max_speed,
                          best$max_gaps, best$quality_pctile * 100),
         description = c("Min detections per track",
                          "Max per-track mean speed (µm/min)",
                          "Max per-track peak speed (µm/min)",
                          "Max frame gaps per track",
                          "Bottom % quality removed")),
  file.path(OUTPUT_DIR, "recommended_params.csv"))

# Ground truth
write_csv(cells_corrected, file.path(OUTPUT_DIR, "ground_truth_cells_corrected.csv"))
write_csv(drift_per_frame, file.path(OUTPUT_DIR, "ground_truth_ysl_drift.csv"))
write_csv(
  tibble(metric = c("inst_speed_median", "inst_speed_mean", "track_mean_speed_median"),
         value  = c(GT$inst_speed_median, GT$inst_speed_mean, GT$track_mean_speed_median),
         unit   = "um_per_min"),
  file.path(OUTPUT_DIR, "ground_truth_summary.csv"))

cat(sprintf("  Exported to %s/:\n", OUTPUT_DIR))
cat(sprintf("    filtered_tracks.csv            — %s tracks\n",
            format(nrow(tracks_out), big.mark = ",")))
cat(sprintf("    filtered_spots.csv             — %s spots\n",
            format(nrow(spots_out), big.mark = ",")))
cat("    sweep_results.csv              — all scored parameter combos\n")
cat("    recommended_params.csv         — chosen filter parameters\n")
cat("    ground_truth_cells_corrected.csv\n")
cat("    ground_truth_ysl_drift.csv\n")
cat("    ground_truth_summary.csv\n")

cat("\n")
cat(strrep("=", 70), "\n")
cat(sprintf("  DONE — %s tracks retained (%.1f%%), score %.3f\n",
            format(best$n_tracks, big.mark = ","),
            best$frac_tracks * 100, best$score))
cat("  Next: source('trackmate_analysis.R')\n")
cat(strrep("=", 70), "\n")
