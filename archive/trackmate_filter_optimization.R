# =============================================================================
# TrackMate Filter Optimization — Finding Params That Match Ground Truth
# =============================================================================
#
# PROBLEM:
#   TrackMate tracking quality depends heavily on QC filter parameters.
#   We have manual ground truth (drift-corrected, 40 cells) giving us target:
#     - Inst. speed median:      1.86 µm/min
#     - Inst. speed mean:        2.34 µm/min
#     - Per-track mean speed:    2.38 µm/min  (median across tracks)
#   TrackMate unfiltered is too high (3.39 µm/min per-track mean).
#   Imaris is closer (2.03) with 5981 tracks.
#
# APPROACH:
#   1. Load TrackMate data once, compute instantaneous speeds per spot
#   2. Sweep over filter parameter combinations
#   3. For each combo, apply filters → compute target metrics
#   4. Score how close each combo is to ground truth
#   5. Show the best-performing parameter sets
#   6. Visualize trade-offs
#
# FILTER PARAMETERS TO SWEEP:
#   - MIN_SPOTS:    minimum number of spots per track (longer = more reliable)
#   - MIN_DURATION: minimum track duration in frames
#   - MAX_MEAN_SPEED: maximum per-track mean speed (µm/min) — removes noise
#   - MAX_MAX_SPEED:  maximum per-track max speed — removes bursty errors
#   - MIN_QUALITY: minimum mean quality score (percentile)
#   - MAX_GAPS:    maximum number of gaps in a track
#   - MIN_CONFINEMENT (optional): remove very confined tracks
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

FRAME_INTERVAL_SEC <- 30
FRAME_INTERVAL_MIN <- FRAME_INTERVAL_SEC / 60
SPEED_CONVERSION   <- 60 / FRAME_INTERVAL_SEC  # µm/frame → µm/min
OUTPUT_DIR <- "analysis_output_filter_sweep"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

theme_sweep <- function(base_size = 10) {
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
cat("  TRACKMATE FILTER OPTIMIZATION\n")
cat(strrep("=", 70), "\n")

# =============================================================================
# GROUND TRUTH TARGETS (from manual_tracking_analysis.R)
# =============================================================================

GT <- list(
  inst_speed_median      = 1.86,
  inst_speed_mean        = 2.34,
  track_mean_speed_median = 2.38
)

cat(sprintf("\n  Ground truth targets:\n"))
cat(sprintf("    Inst. speed median: %.2f µm/min\n", GT$inst_speed_median))
cat(sprintf("    Inst. speed mean:   %.2f µm/min\n", GT$inst_speed_mean))
cat(sprintf("    Per-track speed:    %.2f µm/min\n", GT$track_mean_speed_median))

# =============================================================================
# STEP 1: LOAD TRACKMATE DATA (once)
# =============================================================================

cat("\n=== STEP 1: LOADING TRACKMATE DATA ===\n")

tm_tracks_raw <- read_csv("4D_tracks.csv",
                           show_col_types = FALSE)[-c(1:3), ]
tm_spots_raw  <- read_csv("4D_spots.csv",
                           show_col_types = FALSE)[-c(1:3), ]

tm_tracks <- tm_tracks_raw %>%
  mutate(across(-LABEL, ~as.numeric(.))) %>%
  mutate(
    DURATION_MIN      = TRACK_DURATION * FRAME_INTERVAL_MIN,
    SPEED_MEAN_UM_MIN = TRACK_MEAN_SPEED * SPEED_CONVERSION,
    SPEED_MAX_UM_MIN  = TRACK_MAX_SPEED  * SPEED_CONVERSION,
    SPEED_STD_UM_MIN  = TRACK_STD_SPEED  * SPEED_CONVERSION
  )

tm_spots <- tm_spots_raw %>%
  mutate(across(c(ID, TRACK_ID, QUALITY, POSITION_X, POSITION_Y, POSITION_Z,
                  FRAME), ~as.numeric(.))) %>%
  filter(!is.na(TRACK_ID))

cat(sprintf("  Loaded: %d tracks, %d spots\n", nrow(tm_tracks), nrow(tm_spots)))

# Pre-compute instantaneous speed for ALL spots (this is expensive, do once)
cat("  Computing instantaneous speeds for all spots...\n")
tm_spots_vel <- tm_spots %>%
  arrange(TRACK_ID, FRAME) %>%
  group_by(TRACK_ID) %>%
  mutate(
    dx = POSITION_X - lag(POSITION_X),
    dy = POSITION_Y - lag(POSITION_Y),
    dz = POSITION_Z - lag(POSITION_Z),
    dt_frames = FRAME - lag(FRAME),
    inst_speed = sqrt(dx^2 + dy^2 + dz^2) / (dt_frames * FRAME_INTERVAL_MIN)
  ) %>%
  ungroup() %>%
  filter(!is.na(inst_speed) & dt_frames > 0)

cat(sprintf("  Computed %d instantaneous speed values\n", nrow(tm_spots_vel)))

# =============================================================================
# STEP 2: UNDERSTAND THE DATA BEFORE SWEEPING
# =============================================================================

cat("\n=== STEP 2: DATA OVERVIEW ===\n")

cat(sprintf("\n  Track property ranges:\n"))
cat(sprintf("    NUMBER_SPOTS:   %d – %d  (median: %.0f)\n",
            min(tm_tracks$NUMBER_SPOTS), max(tm_tracks$NUMBER_SPOTS),
            median(tm_tracks$NUMBER_SPOTS)))
cat(sprintf("    DURATION_MIN:   %.1f – %.1f  (median: %.1f)\n",
            min(tm_tracks$DURATION_MIN), max(tm_tracks$DURATION_MIN),
            median(tm_tracks$DURATION_MIN)))
cat(sprintf("    SPEED_MEAN:     %.3f – %.3f  (median: %.3f µm/min)\n",
            min(tm_tracks$SPEED_MEAN_UM_MIN, na.rm=T),
            max(tm_tracks$SPEED_MEAN_UM_MIN, na.rm=T),
            median(tm_tracks$SPEED_MEAN_UM_MIN, na.rm=T)))
cat(sprintf("    SPEED_MAX:      %.3f – %.3f  (median: %.3f µm/min)\n",
            min(tm_tracks$SPEED_MAX_UM_MIN, na.rm=T),
            max(tm_tracks$SPEED_MAX_UM_MIN, na.rm=T),
            median(tm_tracks$SPEED_MAX_UM_MIN, na.rm=T)))
cat(sprintf("    N_GAPS:         %d – %d  (median: %.0f)\n",
            min(tm_tracks$NUMBER_GAPS), max(tm_tracks$NUMBER_GAPS),
            median(tm_tracks$NUMBER_GAPS)))
cat(sprintf("    QUALITY:        %.0f – %.0f  (median: %.0f)\n",
            min(tm_tracks$TRACK_MEAN_QUALITY, na.rm=T),
            max(tm_tracks$TRACK_MEAN_QUALITY, na.rm=T),
            median(tm_tracks$TRACK_MEAN_QUALITY, na.rm=T)))

# =============================================================================
# STEP 3: PARAMETER SWEEP
# =============================================================================

cat("\n=== STEP 3: PARAMETER SWEEP ===\n")

# Define grid
# Key insight: the main problem is short noisy tracks inflating speed.
# Most impactful filters are MIN_SPOTS and MAX_MEAN_SPEED.
sweep_grid <- expand_grid(
  min_spots      = c(3, 5, 10, 15, 20, 30, 50, 75, 100),
  max_mean_speed = c(3, 4, 5, 7, 10, Inf),
  max_max_speed  = c(5, 8, 10, 15, Inf),
  max_gaps       = c(0, 1, 2, 5, Inf),
  quality_pctile = c(0, 0.1, 0.2, 0.3)  # bottom X% removed
)

# Reduce grid size: total = 9 × 6 × 5 × 5 × 4 = 5400 combos
# But the inst. speed computation (filtering spots) is expensive.
# Optimise: pre-index spot velocities by TRACK_ID for fast lookup.
cat(sprintf("  Grid size: %d combinations\n", nrow(sweep_grid)))

# Compute quality thresholds once
quality_vals <- quantile(tm_tracks$TRACK_MEAN_QUALITY,
                          c(0, 0.1, 0.2, 0.3), na.rm = TRUE)
names(quality_vals) <- c("0", "0.1", "0.2", "0.3")

# Pre-split spot velocities by TRACK_ID for fast lookup
cat("  Pre-indexing spot velocities by TRACK_ID...\n")
vel_by_track <- split(tm_spots_vel$inst_speed, tm_spots_vel$TRACK_ID)

# Sweep function
eval_combo <- function(min_spots, max_mean_speed, max_max_speed, max_gaps, quality_pctile) {
  q_thresh <- quality_vals[as.character(quality_pctile)]
  if (is.na(q_thresh)) q_thresh <- 0

  keep <- tm_tracks %>%
    filter(
      NUMBER_SPOTS >= min_spots,
      SPEED_MEAN_UM_MIN <= max_mean_speed,
      SPEED_MAX_UM_MIN  <= max_max_speed,
      NUMBER_GAPS <= max_gaps,
      TRACK_MEAN_QUALITY >= q_thresh
    )

  n_tracks <- nrow(keep)

  if (n_tracks < 10) {
    return(tibble(
      n_tracks = n_tracks, n_spots = 0L,
      inst_speed_median = NA_real_, inst_speed_mean = NA_real_,
      track_speed_median = NA_real_, track_speed_mean = NA_real_,
      duration_median = NA_real_, confinement_median = NA_real_
    ))
  }

  kept_ids <- as.character(keep$TRACK_ID)

  # Track-level metrics
  track_speed_med  <- median(keep$SPEED_MEAN_UM_MIN, na.rm = TRUE)
  track_speed_mean_val <- mean(keep$SPEED_MEAN_UM_MIN, na.rm = TRUE)
  dur_med          <- median(keep$DURATION_MIN, na.rm = TRUE)
  conf_med         <- median(keep$CONFINEMENT_RATIO, na.rm = TRUE)

  # Instantaneous speeds via pre-indexed lookup (much faster than filter)
  inst_speeds <- unlist(vel_by_track[kept_ids], use.names = FALSE)
  n_sp <- length(inst_speeds)
  inst_med  <- if (n_sp > 0) median(inst_speeds, na.rm = TRUE) else NA_real_
  inst_mean <- if (n_sp > 0) mean(inst_speeds, na.rm = TRUE) else NA_real_

  tibble(
    n_tracks = n_tracks, n_spots = n_sp,
    inst_speed_median  = inst_med, inst_speed_mean  = inst_mean,
    track_speed_median = track_speed_med, track_speed_mean = track_speed_mean_val,
    duration_median    = dur_med, confinement_median = conf_med
  )
}

# Run sweep with progress
cat("  Sweeping (this may take a few minutes)...\n")
pb <- txtProgressBar(min = 0, max = nrow(sweep_grid), style = 3)

sweep_list <- vector("list", nrow(sweep_grid))
for (i in seq_len(nrow(sweep_grid))) {
  sweep_list[[i]] <- eval_combo(
    sweep_grid$min_spots[i],
    sweep_grid$max_mean_speed[i],
    sweep_grid$max_max_speed[i],
    sweep_grid$max_gaps[i],
    sweep_grid$quality_pctile[i]
  )
  setTxtProgressBar(pb, i)
}
close(pb)

sweep_metrics <- bind_rows(sweep_list)
sweep_results <- bind_cols(sweep_grid, sweep_metrics)

cat(sprintf("\n  Sweep complete: %d combos evaluated\n", nrow(sweep_results)))
cat(sprintf("  Non-trivial results (≥10 tracks): %d\n",
            sum(!is.na(sweep_results$inst_speed_median))))

# =============================================================================
# STEP 4: SCORE EACH COMBINATION & SELECT BEST
# =============================================================================
#
# SCORING LOGIC:
#   Score = weighted relative error vs ground truth:
#     score = 0.4 × |ISmed - GT| / GT
#           + 0.3 × |ISmean - GT| / GT
#           + 0.3 × |TSmed - GT| / GT
#
# SELECTION METHOD (choose one):
#
#   "pareto"  — Pareto-optimal knee: normalise score & n_tracks to [0,1],
#               find the Pareto-front point closest to the ideal corner
#               (score=0, tracks=max). Balanced trade-off.
#
#   "score"   — Score-priority: keep all combos whose weighted-average
#               relative error vs GT speeds is ≤ SCORE_TOLERANCE (e.g. 5%),
#               then among those pick the one with the most tracks.
#               Strictly prioritises ground truth similarity.
#
# ── CHANGE THIS TO SWITCH ──────────────────────────────────────────────────
SELECTION_METHOD  <- "score"   # "pareto" or "score"
SCORE_TOLERANCE   <- 0.05     # only used when SELECTION_METHOD = "score"
#                                fraction of GT speed: 0.05 = within 5% of GT
# =============================================================================

cat("\n=== STEP 4: SCORING & SELECTION ===\n")
cat(sprintf("  Selection method: %s\n", SELECTION_METHOD))

N_TOTAL_TRACKS <- nrow(tm_tracks)

sweep_scored <- sweep_results %>%
  filter(!is.na(inst_speed_median)) %>%
  mutate(
    # Relative errors (how far each metric is from ground truth)
    err_inst_median = abs(inst_speed_median - GT$inst_speed_median) / GT$inst_speed_median,
    err_inst_mean   = abs(inst_speed_mean   - GT$inst_speed_mean)   / GT$inst_speed_mean,
    err_track_speed = abs(track_speed_median - GT$track_mean_speed_median) / GT$track_mean_speed_median,

    # Combined score (lower = closer to ground truth)
    score = 0.4 * err_inst_median + 0.3 * err_inst_mean + 0.3 * err_track_speed,

    # Fraction of tracks retained (0–1)
    frac_tracks = n_tracks / N_TOTAL_TRACKS
  ) %>%
  arrange(score)

# --- Pareto front (always computed, used for visualisation + optionally for selection) ---
identify_pareto <- function(df) {
  df <- df %>% arrange(score)
  is_pareto <- rep(TRUE, nrow(df))
  max_tracks_so_far <- 0
  for (i in seq_len(nrow(df))) {
    if (df$n_tracks[i] >= max_tracks_so_far) {
      max_tracks_so_far <- df$n_tracks[i]
    } else {
      is_pareto[i] <- FALSE
    }
  }
  df$is_pareto <- is_pareto
  df
}

sweep_scored <- identify_pareto(sweep_scored)
pareto_front <- sweep_scored %>% filter(is_pareto) %>% arrange(score)

cat(sprintf("  Total valid combos: %d\n", nrow(sweep_scored)))
cat(sprintf("  Pareto-optimal combos: %d\n", nrow(pareto_front)))

# --- SELECT BEST COMBO ---
if (SELECTION_METHOD == "pareto") {
  # Pareto knee: normalise both axes to [0,1], find point closest to ideal (score=0, tracks=1)
  pareto_norm <- pareto_front %>%
    mutate(
      score_norm  = (score - min(score)) / (max(score) - min(score) + 1e-9),
      tracks_norm = (n_tracks - min(n_tracks)) / (max(n_tracks) - min(n_tracks) + 1e-9),
      dist_to_ideal = sqrt(score_norm^2 + (1 - tracks_norm)^2)
    )
  best_idx <- which.min(pareto_norm$dist_to_ideal)
  best <- pareto_norm[best_idx, ]
  score_cutoff <- NA  # not used in this mode

  cat(sprintf("\n  Pareto knee selected: score=%.4f with %s tracks\n",
              best$score, format(best$n_tracks, big.mark = ",")))

} else {
  # Score-priority: keep combos within SCORE_TOLERANCE of GT, pick most tracks
  # The score IS the weighted avg relative error vs GT speeds,
  # so score ≤ 0.05 means speeds are on avg within 5% of GT.
  score_cutoff <- SCORE_TOLERANCE

  near_best <- sweep_scored %>%
    filter(score <= score_cutoff) %>%
    arrange(desc(n_tracks))

  if (nrow(near_best) == 0) {
    cat(sprintf("\n  WARNING: No combos within %.0f%% of GT. Relaxing to best available.\n",
                SCORE_TOLERANCE * 100))
    best <- sweep_scored[1, ]  # already sorted by score
    score_cutoff <- best$score
  } else {
    best <- near_best[1, ]
  }

  cat(sprintf("\n  GT speed tolerance:    ±%.0f%% (score ≤ %.4f)\n",
              SCORE_TOLERANCE * 100, score_cutoff))
  cat(sprintf("  Best score found:      %.4f\n", min(sweep_scored$score)))
  cat(sprintf("  Combos within tolerance: %d\n", nrow(near_best)))
  cat(sprintf("  Selected: score=%.4f with %s tracks (most tracks within tolerance)\n",
              best$score, format(best$n_tracks, big.mark = ",")))
}

cat(sprintf("\n  RECOMMENDED PARAMETERS (%s):\n",
            if (SELECTION_METHOD == "pareto") "Pareto knee — balanced trade-off"
            else sprintf("score-priority — within %.0f%% of GT speeds", SCORE_TOLERANCE * 100)))
cat(sprintf("  ┌─────────────────────────────────────────────────────────────┐\n"))
cat(sprintf("  │  min_spots       = %-5d  (min detections per track)       │\n", best$min_spots))
cat(sprintf("  │  max_mean_speed  = %-5.1f  µm/min (removes noisy tracks)  │\n", best$max_mean_speed))
cat(sprintf("  │  max_max_speed   = %-5.1f  µm/min (removes spike jumps)   │\n", best$max_max_speed))
cat(sprintf("  │  max_gaps        = %-5.0f  (max frame gaps per track)      │\n", best$max_gaps))
cat(sprintf("  │  quality_pctile  = %-5.0f%% (bottom quality removed)       │\n", best$quality_pctile * 100))
cat(sprintf("  └─────────────────────────────────────────────────────────────┘\n"))
cat(sprintf("\n  HOW CLOSE TO GROUND TRUTH:\n"))
cat(sprintf("                        Optimized   Ground Truth   Error\n"))
cat(sprintf("    Inst speed median:  %6.3f      %6.3f         %5.1f%%\n",
            best$inst_speed_median, GT$inst_speed_median, best$err_inst_median * 100))
cat(sprintf("    Inst speed mean:    %6.3f      %6.3f         %5.1f%%\n",
            best$inst_speed_mean, GT$inst_speed_mean, best$err_inst_mean * 100))
cat(sprintf("    Track speed median: %6.3f      %6.3f         %5.1f%%\n",
            best$track_speed_median, GT$track_mean_speed_median, best$err_track_speed * 100))
cat(sprintf("\n  TRACKS RETAINED: %d / %d (%.1f%%)\n",
            best$n_tracks, N_TOTAL_TRACKS, best$frac_tracks * 100))

# =============================================================================
# STEP 5: APPLY BEST FILTERS → GET FILTERED DATA
# =============================================================================

cat("\n=== STEP 5: APPLYING BEST FILTERS ===\n")

q_thresh_best <- quantile(tm_tracks$TRACK_MEAN_QUALITY, best$quality_pctile, na.rm = TRUE)

tm_filtered <- tm_tracks %>%
  filter(
    NUMBER_SPOTS >= best$min_spots,
    SPEED_MEAN_UM_MIN <= best$max_mean_speed,
    SPEED_MAX_UM_MIN  <= best$max_max_speed,
    NUMBER_GAPS <= best$max_gaps,
    TRACK_MEAN_QUALITY >= q_thresh_best
  )

tm_filtered_spots <- tm_spots %>%
  filter(TRACK_ID %in% tm_filtered$TRACK_ID)

tm_filtered_vel <- tm_spots_vel %>%
  filter(TRACK_ID %in% tm_filtered$TRACK_ID)

# Also get the default-filtered version for comparison
tm_default <- tm_tracks %>%
  filter(NUMBER_SPOTS >= 5, SPEED_MEAN_UM_MIN <= 10)
tm_default_vel <- tm_spots_vel %>%
  filter(TRACK_ID %in% tm_default$TRACK_ID)

cat(sprintf("  Filtered tracks: %d  (from %d original)\n", nrow(tm_filtered), N_TOTAL_TRACKS))
cat(sprintf("  Filtered spots:  %d  (from %d original)\n", nrow(tm_filtered_spots), nrow(tm_spots)))
cat(sprintf("  Default-filter:  %d tracks (≥5 spots, ≤10 µm/min)\n", nrow(tm_default)))

# =============================================================================
# STEP 6: LOAD GROUND TRUTH & IMARIS FOR COMPARISON
# =============================================================================

cat("\n=== STEP 6: LOADING COMPARISON DATA ===\n")

# --- Manual ground truth (drift-corrected) ---
manual_speeds_path <- "analysis_output_manual/manual_cells_corrected.csv"
if (file.exists(manual_speeds_path)) {
  cat("  Loading pre-computed manual ground truth...\n")
  manual_corr <- read_csv(manual_speeds_path, show_col_types = FALSE)
  manual_speeds <- manual_corr$corr_speed[!is.na(manual_corr$corr_speed)]
} else {
  MANUAL_DIR <- "high_res_mannualtracking"
  if (dir.exists(MANUAL_DIR)) {
    cat("  Loading manual tracking from raw files + YSL drift correction...\n")
    read_manual_pos <- function(folder_name) {
      prefix <- sub("_Statistics$", "", folder_name)
      path   <- file.path(MANUAL_DIR, folder_name, paste0(prefix, "_Position.csv"))
      if (!file.exists(path)) return(NULL)
      df <- read_csv(path, skip = 3, show_col_types = FALSE) %>%
        select(where(~ !all(is.na(.x))))
      if (nrow(df) < 2) return(NULL)
      df %>%
        rename(POSITION_X = `Position X`, POSITION_Y = `Position Y`,
               POSITION_Z = `Position Z`, FRAME = Time) %>%
        mutate(across(c(POSITION_X, POSITION_Y, POSITION_Z, FRAME), as.numeric)) %>%
        select(POSITION_X, POSITION_Y, POSITION_Z, FRAME) %>%
        mutate(cell_name = prefix)
    }
    all_folders <- list.dirs(MANUAL_DIR, recursive = FALSE, full.names = FALSE)
    cell_folders <- grep("^Spots_cell_", all_folders, value = TRUE)
    ysl_folders  <- grep("^Spots_nuclei_ysl_", all_folders, value = TRUE)
    cell_data <- map_dfr(cell_folders, read_manual_pos)
    ysl_data  <- map_dfr(ysl_folders, read_manual_pos)
    ysl_displ <- ysl_data %>%
      arrange(cell_name, FRAME) %>% group_by(cell_name) %>%
      mutate(dx = POSITION_X - lag(POSITION_X),
             dy = POSITION_Y - lag(POSITION_Y),
             dz = POSITION_Z - lag(POSITION_Z)) %>%
      filter(!is.na(dx)) %>% ungroup()
    drift_per_frame <- ysl_displ %>%
      group_by(FRAME) %>%
      summarise(drift_dx = mean(dx), drift_dy = mean(dy), drift_dz = mean(dz),
                .groups = "drop")
    cells_corr <- cell_data %>%
      arrange(cell_name, FRAME) %>% group_by(cell_name) %>%
      mutate(raw_dx = POSITION_X - lag(POSITION_X),
             raw_dy = POSITION_Y - lag(POSITION_Y),
             raw_dz = POSITION_Z - lag(POSITION_Z),
             dt_frames = FRAME - lag(FRAME)) %>%
      filter(!is.na(raw_dx) & dt_frames == 1) %>% ungroup() %>%
      left_join(drift_per_frame, by = "FRAME") %>%
      mutate(corr_dx = raw_dx - drift_dx, corr_dy = raw_dy - drift_dy,
             corr_dz = raw_dz - drift_dz,
             corr_speed = sqrt(corr_dx^2 + corr_dy^2 + corr_dz^2) / FRAME_INTERVAL_MIN)
    manual_speeds <- cells_corr$corr_speed[!is.na(cells_corr$corr_speed)]
    write_csv(cells_corr, file.path(OUTPUT_DIR, "manual_cells_corrected.csv"))
    cat(sprintf("    %d speed values from %d cells\n",
                length(manual_speeds), n_distinct(cells_corr$cell_name)))
  } else {
    cat("  !! Manual tracking folder not found\n")
    manual_speeds <- NULL
  }
}

# --- Imaris ---
imaris_pos_path <- "all_tracks/ome-tiff.companion-bin-8bit-crop_Position.csv"
if (file.exists(imaris_pos_path)) {
  cat("  Loading Imaris data...\n")
  im_raw <- read_csv(imaris_pos_path, skip = 3, show_col_types = FALSE) %>%
    select(where(~ !all(is.na(.x))))
  im_vel <- im_raw %>%
    rename(TRACK_ID = TrackID, POSITION_X = `Position X`,
           POSITION_Y = `Position Y`, POSITION_Z = `Position Z`, FRAME = Time) %>%
    mutate(across(c(TRACK_ID, POSITION_X, POSITION_Y, POSITION_Z, FRAME), as.numeric)) %>%
    arrange(TRACK_ID, FRAME) %>% group_by(TRACK_ID) %>%
    mutate(dx = POSITION_X - lag(POSITION_X), dy = POSITION_Y - lag(POSITION_Y),
           dz = POSITION_Z - lag(POSITION_Z), dt = FRAME - lag(FRAME),
           inst_speed = sqrt(dx^2 + dy^2 + dz^2) / (dt * FRAME_INTERVAL_MIN)) %>%
    ungroup() %>% filter(!is.na(inst_speed) & dt > 0)
  imaris_speeds <- im_vel$inst_speed
  cat(sprintf("    %d speed values from %d tracks\n",
              length(imaris_speeds), n_distinct(im_vel$TRACK_ID)))
} else {
  imaris_speeds <- NULL
}

# =============================================================================
# STEP 7: PLOTS
# =============================================================================

cat("\n=== STEP 7: GENERATING PLOTS ===\n")

all_colors <- c("TM default"   = "grey60",
                "TM optimized" = "#2166AC",
                "Imaris"       = "#FF7F00",
                "Manual (GT)"  = "#4DAF4A")

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 1: Filter Impact — which parameter matters most?
#   One panel per filter parameter, showing how speed changes as you tighten it.
# ═══════════════════════════════════════════════════════════════════════════════

marginal_data <- bind_rows(
  sweep_scored %>%
    group_by(value = min_spots) %>%
    summarise(inst_median = median(inst_speed_median), track_median = median(track_speed_median),
              n_tracks_med = median(n_tracks), .groups = "drop") %>%
    mutate(param = "min_spots"),
  sweep_scored %>% filter(is.finite(max_mean_speed)) %>%
    group_by(value = max_mean_speed) %>%
    summarise(inst_median = median(inst_speed_median), track_median = median(track_speed_median),
              n_tracks_med = median(n_tracks), .groups = "drop") %>%
    mutate(param = "max_mean_speed"),
  sweep_scored %>% filter(is.finite(max_max_speed)) %>%
    group_by(value = max_max_speed) %>%
    summarise(inst_median = median(inst_speed_median), track_median = median(track_speed_median),
              n_tracks_med = median(n_tracks), .groups = "drop") %>%
    mutate(param = "max_max_speed"),
  sweep_scored %>% filter(is.finite(max_gaps)) %>%
    group_by(value = max_gaps) %>%
    summarise(inst_median = median(inst_speed_median), track_median = median(track_speed_median),
              n_tracks_med = median(n_tracks), .groups = "drop") %>%
    mutate(param = "max_gaps"),
  sweep_scored %>%
    group_by(value = quality_pctile * 100) %>%
    summarise(inst_median = median(inst_speed_median), track_median = median(track_speed_median),
              n_tracks_med = median(n_tracks), .groups = "drop") %>%
    mutate(param = "quality_%_removed")
) %>%
  mutate(param = factor(param, levels = c("min_spots", "max_mean_speed", "max_max_speed",
                                           "max_gaps", "quality_%_removed")))

marginal_long <- marginal_data %>%
  pivot_longer(cols = c(inst_median, track_median),
               names_to = "metric", values_to = "speed") %>%
  mutate(metric = recode(metric,
                          "inst_median"  = "Inst. speed median",
                          "track_median" = "Per-track speed median"))

p1_speed <- ggplot(marginal_long, aes(x = value, y = speed, color = metric)) +
  geom_line(linewidth = 1) + geom_point(size = 1.5) +
  geom_hline(yintercept = GT$inst_speed_median, linetype = "dashed", color = "#4DAF4A", linewidth = 0.5) +
  geom_hline(yintercept = GT$track_mean_speed_median, linetype = "dotted", color = "#4DAF4A", linewidth = 0.5) +
  scale_color_manual(values = c("Inst. speed median" = "#2166AC", "Per-track speed median" = "#B2182B")) +
  facet_wrap(~param, scales = "free_x", nrow = 1) +
  labs(title = "How each filter parameter affects speed metrics",
       subtitle = "Green dashed/dotted = ground truth targets | Each panel varies one parameter, averaging over the others",
       x = "Parameter value", y = "Speed (µm/min)", color = NULL) +
  theme_sweep(base_size = 9) +
  theme(strip.text = element_text(size = 8))

p1_n <- ggplot(marginal_data, aes(x = value, y = n_tracks_med)) +
  geom_line(linewidth = 0.8, color = "grey40") + geom_point(size = 1.5) +
  scale_y_log10() +
  facet_wrap(~param, scales = "free_x", nrow = 1) +
  labs(title = "Tracks remaining after each filter",
       x = "Parameter value", y = "N tracks (log scale)") +
  theme_sweep(base_size = 9) +
  theme(strip.text = element_text(size = 8))

p1_combined <- p1_speed / p1_n +
  plot_annotation(
    title = "Filter Impact Analysis",
    subtitle = "Which QC parameters have the biggest effect on speed metrics and data retention?",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50"))
  )

ggsave(file.path(OUTPUT_DIR, "01_filter_impact.pdf"), p1_combined, width = 16, height = 8)
cat("  Saved: 01_filter_impact.pdf\n")

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 2: 2D Heatmaps (MIN_SPOTS × MAX_MEAN_SPEED)
#   Three separate plots (each with own color scale) combined with patchwork.
# ═══════════════════════════════════════════════════════════════════════════════

heatmap_data <- sweep_scored %>%
  filter(is.finite(max_mean_speed)) %>%
  group_by(min_spots, max_mean_speed) %>%
  summarise(mean_score   = mean(score, na.rm = TRUE),
            mean_inst    = mean(inst_speed_median, na.rm = TRUE),
            mean_tracks  = mean(n_tracks, na.rm = TRUE),
            .groups = "drop")

p2_score <- ggplot(heatmap_data, aes(x = factor(min_spots), y = factor(max_mean_speed),
                                      fill = mean_score)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", mean_score)), size = 2.3, color = "white") +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Score") +
  labs(title = "Score (lower = closer to GT)",
       x = "Min spots", y = "Max mean speed (µm/min)") +
  theme_sweep(base_size = 9)

p2_inst <- ggplot(heatmap_data, aes(x = factor(min_spots), y = factor(max_mean_speed),
                                     fill = mean_inst)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", mean_inst)), size = 2.3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = GT$inst_speed_median,
                        name = "µm/min") +
  labs(title = sprintf("Inst. Speed Median (target: %.2f)", GT$inst_speed_median),
       x = "Min spots", y = "Max mean speed (µm/min)") +
  theme_sweep(base_size = 9)

p2_n <- ggplot(heatmap_data, aes(x = factor(min_spots), y = factor(max_mean_speed),
                                  fill = mean_tracks)) +
  geom_tile() +
  geom_text(aes(label = format(round(mean_tracks), big.mark = ",")), size = 2.1) +
  scale_fill_viridis_c(option = "viridis", name = "N tracks",
                        labels = scales::comma) +
  labs(title = "Tracks Retained",
       x = "Min spots", y = "Max mean speed (µm/min)") +
  theme_sweep(base_size = 9)

p2_combined <- p2_score + p2_inst + p2_n +
  plot_annotation(
    title = "Parameter Space: MIN_SPOTS × MAX_MEAN_SPEED",
    subtitle = "Each heatmap has its own color scale | Values averaged over other filter params",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey50"))
  )

ggsave(file.path(OUTPUT_DIR, "02_parameter_heatmaps.pdf"), p2_combined, width = 18, height = 6)
cat("  Saved: 02_parameter_heatmaps.pdf\n")

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 3: Pareto Front — Score vs N Tracks
#   Shows the trade-off and where the recommended combo sits.
# ═══════════════════════════════════════════════════════════════════════════════

p3 <- ggplot(sweep_scored, aes(x = n_tracks, y = score)) +
  geom_point(alpha = 0.08, size = 0.8, color = "grey50") +
  geom_point(data = pareto_front, color = "#E41A1C", size = 2) +
  geom_line(data = pareto_front, color = "#E41A1C", linewidth = 0.8) +
  geom_point(data = best, aes(x = n_tracks, y = score),
             color = "#2166AC", size = 5, shape = 18) +
  {if (!is.na(score_cutoff))
    geom_hline(yintercept = score_cutoff, linetype = "dotted", color = "#2166AC",
               linewidth = 0.5)} +
  {if (!is.na(score_cutoff))
    annotate("text", x = min(sweep_scored$n_tracks), y = score_cutoff + 0.005,
             label = sprintf("within %.0f%% of GT speeds", SCORE_TOLERANCE * 100),
             hjust = 0, size = 2.5, color = "#2166AC")} +
  annotate("label", x = best$n_tracks, y = best$score + 0.02,
           label = sprintf("SELECTED\n%s tracks, score=%.3f",
                           format(best$n_tracks, big.mark = ","), best$score),
           size = 3, color = "#2166AC", fill = "white", label.size = 0.3) +
  scale_x_log10(labels = scales::comma) +
  labs(title = "Pareto Front: Ground Truth Similarity vs Data Retention",
       subtitle = paste0("Red line = Pareto front | ",
                         if (SELECTION_METHOD == "score")
                           sprintf("Blue dotted = within %.0f%% of GT speeds\n", SCORE_TOLERANCE * 100)
                         else "",
                         "Blue diamond = selected (",
                         if (SELECTION_METHOD == "pareto") "Pareto knee"
                         else "within GT tolerance, most tracks", ")"),
       x = "Number of tracks retained (log scale)",
       y = "Score (lower = closer to ground truth)") +
  theme_sweep()

ggsave(file.path(OUTPUT_DIR, "03_pareto_front.pdf"), p3, width = 10, height = 7)
cat("  Saved: 03_pareto_front.pdf\n")

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 4: FACETED speed distributions — one panel per dataset
#   Much cleaner than overlapping 4 densities on one plot.
# ═══════════════════════════════════════════════════════════════════════════════

speed_all <- bind_rows(
  tibble(speed = tm_default_vel$inst_speed, dataset = "TM default"),
  tibble(speed = tm_filtered_vel$inst_speed, dataset = "TM optimized"),
  if (!is.null(imaris_speeds)) tibble(speed = imaris_speeds, dataset = "Imaris") else NULL,
  if (!is.null(manual_speeds)) tibble(speed = manual_speeds, dataset = "Manual (GT)") else NULL
) %>%
  mutate(dataset = factor(dataset, levels = c("TM default", "TM optimized", "Imaris", "Manual (GT)")))

speed_stats <- speed_all %>%
  group_by(dataset) %>%
  summarise(median = median(speed, na.rm = TRUE),
            mean   = mean(speed, na.rm = TRUE),
            n      = n(), .groups = "drop")

p4_facet <- ggplot(speed_all, aes(x = speed, fill = dataset)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, alpha = 0.7) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  geom_vline(data = speed_stats, aes(xintercept = median),
             linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_vline(xintercept = GT$inst_speed_median,
             linetype = "dotted", color = "#4DAF4A", linewidth = 0.6) +
  geom_text(data = speed_stats,
            aes(x = Inf, y = Inf,
                label = sprintf("median = %.2f\nmean = %.2f\nn = %s",
                                median, mean, format(n, big.mark = ","))),
            hjust = 1.1, vjust = 1.3, size = 2.8, inherit.aes = FALSE) +
  facet_wrap(~dataset, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = all_colors) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(title = "Instantaneous Speed Distributions — Each Pipeline Separately",
       subtitle = sprintf("Black dashed = per-dataset median | Green dotted = ground truth target (%.2f µm/min)",
                          GT$inst_speed_median),
       x = "Instantaneous speed (µm/min)", y = "Density") +
  theme_sweep(base_size = 9) +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "04_faceted_speed_distributions.pdf"), p4_facet, width = 16, height = 5)
cat("  Saved: 04_faceted_speed_distributions.pdf\n")

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 5: Overlay comparison — density + box side by side
# ═══════════════════════════════════════════════════════════════════════════════

p5_dens <- ggplot(speed_all, aes(x = speed, color = dataset)) +
  geom_density(linewidth = 0.9) +
  geom_vline(xintercept = GT$inst_speed_median, linetype = "dashed",
             color = "#4DAF4A", linewidth = 0.5) +
  scale_color_manual(values = all_colors) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(title = "Overlay: Density Curves",
       subtitle = "Green dashed = ground truth median",
       x = "Speed (µm/min)", y = "Density", color = NULL) +
  theme_sweep()

p5_box <- ggplot(speed_all, aes(x = dataset, y = speed, fill = dataset)) +
  geom_boxplot(alpha = 0.6, outlier.alpha = 0.02, outlier.size = 0.3) +
  geom_hline(yintercept = GT$inst_speed_median, linetype = "dashed",
             color = "#4DAF4A", linewidth = 0.6) +
  scale_fill_manual(values = all_colors) +
  coord_cartesian(ylim = c(0, 10)) +
  labs(title = "Box Plot Comparison",
       subtitle = "Green dashed = ground truth median",
       x = NULL, y = "Speed (µm/min)") +
  theme_sweep() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 1))

p5 <- p5_dens + p5_box +
  plot_annotation(
    title = "Speed Comparison Across All Pipelines",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "05_overlay_comparison.pdf"), p5, width = 14, height = 6)
cat("  Saved: 05_overlay_comparison.pdf\n")

# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 6: Track duration & quality distributions (before vs after filtering)
# ═══════════════════════════════════════════════════════════════════════════════

dur_compare <- bind_rows(
  tibble(duration = tm_tracks$DURATION_MIN, group = "All tracks (unfiltered)"),
  tibble(duration = tm_default$DURATION_MIN, group = "TM default"),
  tibble(duration = tm_filtered$DURATION_MIN, group = "TM optimized")
) %>%
  mutate(group = factor(group, levels = c("All tracks (unfiltered)", "TM default", "TM optimized")))

p6_dur <- ggplot(dur_compare, aes(x = duration, fill = group)) +
  geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
  facet_wrap(~group, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("All tracks (unfiltered)" = "grey80",
                                "TM default" = "grey50",
                                "TM optimized" = "#2166AC")) +
  labs(title = "Track Duration Distribution — Before vs After Filtering",
       subtitle = "Short noisy tracks are the main source of speed inflation",
       x = "Track duration (min)", y = "Count") +
  theme_sweep(base_size = 9) +
  theme(legend.position = "none")

spots_compare <- bind_rows(
  tibble(n_spots = tm_tracks$NUMBER_SPOTS, group = "All tracks"),
  tibble(n_spots = tm_default$NUMBER_SPOTS, group = "TM default"),
  tibble(n_spots = tm_filtered$NUMBER_SPOTS, group = "TM optimized")
) %>%
  mutate(group = factor(group, levels = c("All tracks", "TM default", "TM optimized")))

p6_spots <- ggplot(spots_compare, aes(x = n_spots, fill = group)) +
  geom_histogram(bins = 60, alpha = 0.7) +
  facet_wrap(~group, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("All tracks" = "grey80", "TM default" = "grey50",
                                "TM optimized" = "#2166AC")) +
  labs(title = "Number of Spots per Track",
       x = "N spots", y = "Count") +
  theme_sweep(base_size = 9) +
  theme(legend.position = "none")

p6 <- p6_dur / p6_spots +
  plot_annotation(
    title = "What the Filters Remove: Duration & Track Length",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "06_filtering_effect.pdf"), p6, width = 16, height = 8)
cat("  Saved: 06_filtering_effect.pdf\n")

# =============================================================================
# STEP 8: CONSOLE SUMMARY
# =============================================================================

cat("\n=== STEP 8: FINAL SUMMARY ===\n\n")

cat("  Speed comparison across pipelines:\n")
cat(sprintf("  %-18s  %8s  %8s  %10s\n", "Pipeline", "Median", "Mean", "N values"))
cat("  ", strrep("-", 50), "\n")
for (i in 1:nrow(speed_stats)) {
  r <- speed_stats[i, ]
  cat(sprintf("  %-18s  %8.3f  %8.3f  %10s\n",
              as.character(r$dataset), r$median, r$mean, format(r$n, big.mark = ",")))
}
cat(sprintf("  %-18s  %8.3f  %8.3f  %10s\n",
            "GROUND TRUTH", GT$inst_speed_median, GT$inst_speed_mean, "—"))

# Top 5 Pareto-optimal alternatives
cat("\n  Top 5 Pareto-optimal alternatives (diverse trade-offs):\n")
top5 <- head(pareto_front, 5)
cat(sprintf("  %-4s  %5s %6s %6s %5s %5s  %6s %6s %6s  %6s  %6s\n",
            "#", "MinSp", "MaxMS", "MaxXS", "MaxG", "Q%",
            "ISmed", "ISmn", "TSmed", "Score", "Ntrk"))
cat("  ", strrep("-", 80), "\n")
for (i in 1:nrow(top5)) {
  r <- top5[i, ]
  marker <- if (i == best_idx && best_idx <= 5) " <-BEST" else ""
  cat(sprintf("  %-4d  %5d %6.1f %6.1f %5.0f %5.0f  %6.3f %6.3f %6.3f  %6.3f  %6d%s\n",
              i, r$min_spots, r$max_mean_speed, r$max_max_speed, r$max_gaps,
              r$quality_pctile * 100,
              r$inst_speed_median, r$inst_speed_mean, r$track_speed_median,
              r$score, r$n_tracks, marker))
}

# =============================================================================
# STEP 9: EXPORT FILTERED DATA
# =============================================================================

cat("\n=== STEP 9: EXPORTING FILTERED DATA ===\n")

# --- 9A: Sweep results ---
write_csv(sweep_scored, file.path(OUTPUT_DIR, "sweep_all_results.csv"))
write_csv(pareto_front, file.path(OUTPUT_DIR, "sweep_pareto_front.csv"))

# --- 9B: Recommended parameters ---
best_params <- tibble(
  parameter   = c("min_spots", "max_mean_speed_um_min", "max_max_speed_um_min",
                   "max_gaps", "quality_bottom_pct_removed"),
  value       = c(best$min_spots, best$max_mean_speed, best$max_max_speed,
                   best$max_gaps, best$quality_pctile * 100),
  description = c("Minimum number of detections (spots) per track",
                   "Maximum per-track mean speed in µm/min — removes erroneously fast tracks",
                   "Maximum per-track peak speed in µm/min — removes single-step spike jumps",
                   "Maximum number of detection gaps allowed per track",
                   "Bottom percentage of quality scores removed (0 = keep all)")
)
write_csv(best_params, file.path(OUTPUT_DIR, "recommended_params.csv"))

# --- 9C: FILTERED TRACKS CSV ---
#   Re-export the original TrackMate tracks CSV format, keeping only passing tracks.
#   We include the original columns plus the computed µm/min speeds.
tracks_out <- tm_filtered %>%
  select(TRACK_ID, LABEL, NUMBER_SPOTS, NUMBER_GAPS, NUMBER_SPLITS, NUMBER_MERGES,
         LONGEST_GAP, TRACK_DURATION, TRACK_START, TRACK_STOP, TRACK_DISPLACEMENT,
         TRACK_X_LOCATION, TRACK_Y_LOCATION, TRACK_Z_LOCATION,
         TRACK_MEAN_SPEED, TRACK_MAX_SPEED, TRACK_MIN_SPEED,
         TRACK_MEDIAN_SPEED, TRACK_STD_SPEED, TRACK_MEAN_QUALITY,
         TOTAL_DISTANCE_TRAVELED, MAX_DISTANCE_TRAVELED,
         CONFINEMENT_RATIO, MEAN_STRAIGHT_LINE_SPEED,
         LINEARITY_OF_FORWARD_PROGRESSION, MEAN_DIRECTIONAL_CHANGE_RATE,
         DURATION_MIN, SPEED_MEAN_UM_MIN, SPEED_MAX_UM_MIN, SPEED_STD_UM_MIN)

write_csv(tracks_out, file.path(OUTPUT_DIR, "filtered_tracks.csv"))

# --- 9D: FILTERED SPOTS CSV ---
#   All spots belonging to the filtered tracks, with original + computed columns.
spots_out <- tm_filtered_spots %>%
  select(any_of(c("ID", "LABEL", "TRACK_ID", "QUALITY",
                   "POSITION_X", "POSITION_Y", "POSITION_Z",
                   "FRAME", "RADIUS", "VISIBILITY",
                   "MANUAL_SPOT_COLOR", "MEAN_INTENSITY_CH1",
                   "MEDIAN_INTENSITY_CH1", "MIN_INTENSITY_CH1",
                   "MAX_INTENSITY_CH1", "TOTAL_INTENSITY_CH1",
                   "STD_INTENSITY_CH1", "CONTRAST_CH1", "SNR_CH1")))

write_csv(spots_out, file.path(OUTPUT_DIR, "filtered_spots.csv"))

cat(sprintf("  Exported to %s/:\n", OUTPUT_DIR))
cat(sprintf("    filtered_tracks.csv       — %s tracks (ready for downstream analysis)\n",
            format(nrow(tracks_out), big.mark = ",")))
cat(sprintf("    filtered_spots.csv        — %s spots from those tracks\n",
            format(nrow(spots_out), big.mark = ",")))
cat("    sweep_all_results.csv     — all %s scored parameter combinations\n")
cat("    sweep_pareto_front.csv    — Pareto-optimal combos only\n")
cat("    recommended_params.csv    — the chosen filter parameters\n")

cat("\n")
cat(strrep("=", 70), "\n")
cat("  FILTER OPTIMIZATION COMPLETE\n")
cat(sprintf("  Recommended: %d tracks retained (%.1f%%) with score %.3f\n",
            best$n_tracks, best$frac_tracks * 100, best$score))
cat(strrep("=", 70), "\n")
