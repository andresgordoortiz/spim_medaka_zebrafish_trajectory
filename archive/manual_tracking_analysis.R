# =============================================================================
# Manual Tracking Analysis — YSL-Corrected Ground Truth
# + 3-Way Comparison with TrackMate & Imaris
# =============================================================================
#
# CONTEXT:
#   - 40 cells + 5 YSL nuclei manually tracked in Imaris (different experiment,
#     same species & conditions → speeds should be comparable)
#   - Frame interval = 30 sec (uncalibrated in Imaris: 1 frame = 1 "second")
#   - The embryo moved substantially during acquisition
#   - YSL nuclei are embedded in the yolk syncytial layer and do NOT migrate
#     → their apparent movement = embryo drift
#   - We subtract the mean YSL displacement at each frame gap to correct all
#     cell tracks for embryo movement
#
# PIPELINE:
#   STEP 1: Load all manual tracks (cells + YSL)
#   STEP 2: Inspect raw data & embryo drift
#   STEP 3: Compute YSL-based drift correction per frame
#   STEP 4: Apply correction, compute drift-corrected velocities
#   STEP 5: Visualize corrected manual tracking
#   STEP 6: Load TrackMate & Imaris summary stats (from previous scripts)
#   STEP 7: 3-way comparison
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

# =============================================================================
# PARAMETERS
# =============================================================================

MANUAL_DIR    <- "high_res_mannualtracking"
FRAME_INTERVAL_SEC <- 30
FRAME_INTERVAL_MIN <- FRAME_INTERVAL_SEC / 60
OUTPUT_DIR    <- "analysis_output_manual"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

software_colors <- c("Manual (corrected)" = "#4DAF4A",
                      "TrackMate"          = "#2166AC",
                      "Imaris"             = "#B2182B")

theme_manual <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 1, color = "grey50"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = base_size - 1),
      strip.text = element_text(face = "bold")
    )
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("  MANUAL TRACKING — YSL DRIFT CORRECTION & 3-WAY COMPARISON\n")
cat(strrep("=", 70), "\n")
cat(sprintf("  Frame interval: %d sec (%.2f min)\n", FRAME_INTERVAL_SEC, FRAME_INTERVAL_MIN))

# #############################################################################
# STEP 1: LOAD ALL MANUAL TRACKS
# #############################################################################

cat("\n=== STEP 1: LOADING MANUAL TRACKS ===\n")

# Helper: read one Imaris Position CSV (skip 3 header lines, drop trailing NA col)
read_manual_position <- function(folder_name) {
  # Derive the file prefix from folder name (strip "_Statistics")
  prefix <- sub("_Statistics$", "", folder_name)
  path   <- file.path(MANUAL_DIR, folder_name, paste0(prefix, "_Position.csv"))

  if (!file.exists(path)) {
    warning(sprintf("Missing: %s", path))
    return(NULL)
  }

  df <- read_csv(path, skip = 3, show_col_types = FALSE) %>%
    select(where(~ !all(is.na(.x))))

  if (nrow(df) < 2) return(NULL)

  df %>%
    rename(POSITION_X = `Position X`,
           POSITION_Y = `Position Y`,
           POSITION_Z = `Position Z`,
           FRAME      = Time) %>%
    mutate(across(c(POSITION_X, POSITION_Y, POSITION_Z, FRAME), as.numeric)) %>%
    select(POSITION_X, POSITION_Y, POSITION_Z, FRAME) %>%
    mutate(cell_name = prefix)
}

# --- Discover folders ---
all_folders <- list.dirs(MANUAL_DIR, recursive = FALSE, full.names = FALSE)
cell_folders <- grep("^Spots_cell_", all_folders, value = TRUE)
ysl_folders  <- grep("^Spots_nuclei_ysl_", all_folders, value = TRUE)

cat(sprintf("  Cell folders: %d\n", length(cell_folders)))
cat(sprintf("  YSL folders:  %d\n", length(ysl_folders)))

# --- Load all ---
cell_data <- map_dfr(cell_folders, read_manual_position) %>%
  mutate(cell_type = "cell")

ysl_data <- map_dfr(ysl_folders, read_manual_position) %>%
  mutate(cell_type = "ysl")

manual_all <- bind_rows(cell_data, ysl_data)

cat(sprintf("  Loaded %d cell tracks (%d timepoints)\n",
            n_distinct(cell_data$cell_name), nrow(cell_data)))
cat(sprintf("  Loaded %d YSL tracks  (%d timepoints)\n",
            n_distinct(ysl_data$cell_name), nrow(ysl_data)))

# Track lengths
track_lengths <- manual_all %>%
  group_by(cell_name, cell_type) %>%
  summarise(n_frames = n(),
            first_frame = min(FRAME),
            last_frame  = max(FRAME), .groups = "drop") %>%
  arrange(cell_type, desc(n_frames))

cat("\n  Track lengths summary:\n")
cat(sprintf("    Cells — median: %.0f frames | range: %d–%d\n",
            median(track_lengths$n_frames[track_lengths$cell_type == "cell"]),
            min(track_lengths$n_frames[track_lengths$cell_type == "cell"]),
            max(track_lengths$n_frames[track_lengths$cell_type == "cell"])))
cat(sprintf("    YSL   — median: %.0f frames | range: %d–%d\n",
            median(track_lengths$n_frames[track_lengths$cell_type == "ysl"]),
            min(track_lengths$n_frames[track_lengths$cell_type == "ysl"]),
            max(track_lengths$n_frames[track_lengths$cell_type == "ysl"])))

# #############################################################################
# STEP 2: INSPECT RAW DATA & EMBRYO DRIFT
# #############################################################################

cat("\n=== STEP 2: RAW DATA & EMBRYO DRIFT ===\n")

# Raw XY tracks — all cells + YSL
p_raw_xy <- ggplot(manual_all, aes(x = POSITION_X, y = POSITION_Y,
                                    group = cell_name, color = cell_type)) +
  geom_path(linewidth = 0.5, alpha = 0.7) +
  geom_point(data = manual_all %>% group_by(cell_name) %>% slice(1),
             size = 2, shape = 16) +
  scale_color_manual(values = c("cell" = "#2166AC", "ysl" = "#E66101"),
                     labels = c("Cell", "YSL (reference)")) +
  coord_fixed() +
  labs(title = "Raw XY Tracks (before drift correction)",
       subtitle = "YSL nuclei should not move — apparent motion = embryo drift",
       x = "X (µm)", y = "Y (µm)", color = NULL) +
  theme_manual()

# YSL trajectories alone — to visualize the drift
p_ysl_drift <- ggplot(ysl_data, aes(x = POSITION_X, y = POSITION_Y,
                                     group = cell_name, color = cell_name)) +
  geom_path(linewidth = 0.8, alpha = 0.8) +
  geom_point(data = ysl_data %>% group_by(cell_name) %>% slice(1),
             size = 3, shape = 17) +
  labs(title = "YSL Nuclei Trajectories",
       subtitle = "All apparent motion = embryo drift",
       x = "X (µm)", y = "Y (µm)", color = NULL) +
  theme_manual()

# Per-frame YSL displacement (raw, before averaging)
ysl_displacements <- ysl_data %>%
  arrange(cell_name, FRAME) %>%
  group_by(cell_name) %>%
  mutate(
    dx = POSITION_X - lag(POSITION_X),
    dy = POSITION_Y - lag(POSITION_Y),
    dz = POSITION_Z - lag(POSITION_Z),
    displacement = sqrt(dx^2 + dy^2 + dz^2),
    time_min = FRAME * FRAME_INTERVAL_MIN
  ) %>%
  ungroup() %>%
  filter(!is.na(dx))

p_ysl_disp_time <- ggplot(ysl_displacements,
                           aes(x = time_min, y = displacement, color = cell_name)) +
  geom_line(alpha = 0.8) +
  labs(title = "YSL Per-Step Displacement Over Time",
       subtitle = "Should be ~constant if drift is uniform",
       x = "Time (min)", y = "Step displacement (µm)", color = NULL) +
  theme_manual()

# Per-axis drift over time for YSL
ysl_pos_long <- ysl_data %>%
  mutate(time_min = FRAME * FRAME_INTERVAL_MIN) %>%
  pivot_longer(cols = c(POSITION_X, POSITION_Y, POSITION_Z),
               names_to = "axis", values_to = "position") %>%
  mutate(axis = sub("POSITION_", "", axis))

p_ysl_axes <- ggplot(ysl_pos_long, aes(x = time_min, y = position,
                                        group = cell_name, color = cell_name)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~axis, scales = "free_y", ncol = 1) +
  labs(title = "YSL Position Over Time (per axis)",
       subtitle = "Coherent trends = embryo drift",
       x = "Time (min)", y = "Position (µm)", color = NULL) +
  theme_manual()

raw_overview <- (p_raw_xy + p_ysl_drift) / (p_ysl_disp_time + p_ysl_axes) +
  plot_annotation(
    title = "Step 2: Raw Data & Embryo Drift Inspection",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "01_raw_data_drift.pdf"), raw_overview,
       width = 16, height = 14)
cat("  Saved: 01_raw_data_drift.pdf\n")

# #############################################################################
# STEP 3: COMPUTE YSL-BASED DRIFT CORRECTION
# #############################################################################

cat("\n=== STEP 3: YSL DRIFT CORRECTION ===\n")

# Strategy: At each frame gap (f → f+1), compute the mean displacement
# across all YSL nuclei. This is the "embryo movement" for that interval.
# We subtract it from every cell's displacement.

# Per-frame per-YSL displacements
ysl_frame_disp <- ysl_data %>%
  arrange(cell_name, FRAME) %>%
  group_by(cell_name) %>%
  mutate(
    dx = POSITION_X - lag(POSITION_X),
    dy = POSITION_Y - lag(POSITION_Y),
    dz = POSITION_Z - lag(POSITION_Z)
  ) %>%
  ungroup() %>%
  filter(!is.na(dx))

# Average drift per frame across YSL nuclei
drift_per_frame <- ysl_frame_disp %>%
  group_by(FRAME) %>%
  summarise(
    drift_dx = mean(dx, na.rm = TRUE),
    drift_dy = mean(dy, na.rm = TRUE),
    drift_dz = mean(dz, na.rm = TRUE),
    drift_sd_x = sd(dx, na.rm = TRUE),
    drift_sd_y = sd(dy, na.rm = TRUE),
    drift_sd_z = sd(dz, na.rm = TRUE),
    n_ysl = n(),
    .groups = "drop"
  ) %>%
  mutate(
    drift_magnitude = sqrt(drift_dx^2 + drift_dy^2 + drift_dz^2),
    time_min = FRAME * FRAME_INTERVAL_MIN
  )

cat(sprintf("  Drift correction computed for %d frame gaps\n", nrow(drift_per_frame)))
cat(sprintf("  Mean drift magnitude: %.3f µm/frame (%.3f µm/min)\n",
            mean(drift_per_frame$drift_magnitude),
            mean(drift_per_frame$drift_magnitude) / FRAME_INTERVAL_MIN))
cat(sprintf("  Max drift magnitude:  %.3f µm/frame\n",
            max(drift_per_frame$drift_magnitude)))
cat(sprintf("  YSL nuclei per frame: %d–%d\n",
            min(drift_per_frame$n_ysl), max(drift_per_frame$n_ysl)))

# Visualize drift vector over time
p_drift_mag <- ggplot(drift_per_frame, aes(x = time_min, y = drift_magnitude)) +
  geom_line(color = "#E66101", linewidth = 0.8) +
  geom_point(color = "#E66101", size = 1) +
  labs(title = "Embryo Drift Magnitude Per Frame",
       subtitle = "Computed from mean YSL displacement",
       x = "Time (min)", y = "Drift (µm/frame)") +
  theme_manual()

drift_long <- drift_per_frame %>%
  pivot_longer(cols = c(drift_dx, drift_dy, drift_dz),
               names_to = "axis", values_to = "drift") %>%
  mutate(axis = sub("drift_d", "", axis) %>% toupper())

p_drift_axes <- ggplot(drift_long, aes(x = time_min, y = drift, color = axis)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("X" = "#E41A1C", "Y" = "#377EB8", "Z" = "#4DAF4A")) +
  labs(title = "Embryo Drift Per Axis",
       subtitle = "Mean YSL displacement per frame gap",
       x = "Time (min)", y = "Drift (µm)", color = "Axis") +
  theme_manual()

# Cumulative drift path
cum_drift <- drift_per_frame %>%
  arrange(FRAME) %>%
  mutate(
    cum_x = cumsum(drift_dx),
    cum_y = cumsum(drift_dy),
    cum_z = cumsum(drift_dz)
  )

p_drift_path <- ggplot(cum_drift, aes(x = cum_x, y = cum_y, color = time_min)) +
  geom_path(linewidth = 1.2) +
  geom_point(size = 1) +
  scale_color_viridis_c(name = "Time\n(min)") +
  coord_fixed() +
  labs(title = "Cumulative Embryo Drift (XY)",
       subtitle = "Total path the embryo drifted",
       x = "Cumulative X drift (µm)", y = "Cumulative Y drift (µm)") +
  theme_manual()

drift_combined <- (p_drift_mag + p_drift_axes) / (p_drift_path + plot_spacer()) +
  plot_annotation(
    title = "Step 3: YSL-Based Drift Correction",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "02_drift_correction.pdf"), drift_combined,
       width = 14, height = 10)
cat("  Saved: 02_drift_correction.pdf\n")

# #############################################################################
# STEP 4: APPLY CORRECTION & COMPUTE VELOCITIES
# #############################################################################

cat("\n=== STEP 4: DRIFT-CORRECTED VELOCITIES ===\n")

# For each cell, compute raw displacements and subtract drift
correct_cell <- function(df, drift_df) {
  df %>%
    arrange(FRAME) %>%
    mutate(
      # Raw displacements
      raw_dx = POSITION_X - lag(POSITION_X),
      raw_dy = POSITION_Y - lag(POSITION_Y),
      raw_dz = POSITION_Z - lag(POSITION_Z),
      dt_frames = FRAME - lag(FRAME)
    ) %>%
    filter(!is.na(raw_dx) & dt_frames == 1) %>%
    left_join(drift_df %>% select(FRAME, drift_dx, drift_dy, drift_dz),
              by = "FRAME") %>%
    mutate(
      # Corrected displacements = raw - drift
      corr_dx = raw_dx - drift_dx,
      corr_dy = raw_dy - drift_dy,
      corr_dz = raw_dz - drift_dz,

      # Speeds (µm/min)
      raw_speed  = sqrt(raw_dx^2  + raw_dy^2  + raw_dz^2)  / FRAME_INTERVAL_MIN,
      corr_speed = sqrt(corr_dx^2 + corr_dy^2 + corr_dz^2) / FRAME_INTERVAL_MIN,
      drift_speed = sqrt(drift_dx^2 + drift_dy^2 + drift_dz^2) / FRAME_INTERVAL_MIN,

      time_min = FRAME * FRAME_INTERVAL_MIN
    )
}

cells_corrected <- cell_data %>%
  group_by(cell_name) %>%
  group_modify(~ correct_cell(.x, drift_per_frame)) %>%
  ungroup()

# Also correct YSL to verify (corrected YSL should be ~0 speed)
ysl_corrected <- ysl_data %>%
  group_by(cell_name) %>%
  group_modify(~ correct_cell(.x, drift_per_frame)) %>%
  ungroup()

# Summary
cat("\n  Speed Statistics (µm/min):\n")
cat(sprintf("    %-25s  %8s  %8s  %8s\n", "", "Median", "Mean", "SD"))

raw_med  <- median(cells_corrected$raw_speed, na.rm = TRUE)
raw_mean <- mean(cells_corrected$raw_speed, na.rm = TRUE)
raw_sd   <- sd(cells_corrected$raw_speed, na.rm = TRUE)
cat(sprintf("    %-25s  %8.3f  %8.3f  %8.3f\n", "Cells (raw)", raw_med, raw_mean, raw_sd))

corr_med  <- median(cells_corrected$corr_speed, na.rm = TRUE)
corr_mean <- mean(cells_corrected$corr_speed, na.rm = TRUE)
corr_sd   <- sd(cells_corrected$corr_speed, na.rm = TRUE)
cat(sprintf("    %-25s  %8.3f  %8.3f  %8.3f\n", "Cells (corrected)", corr_med, corr_mean, corr_sd))

drift_med <- median(cells_corrected$drift_speed, na.rm = TRUE)
cat(sprintf("    %-25s  %8.3f\n", "Drift speed (median)", drift_med))

ysl_corr_med <- median(ysl_corrected$corr_speed, na.rm = TRUE)
cat(sprintf("    %-25s  %8.3f  (should be ~0)\n", "YSL (corrected)", ysl_corr_med))

# Per-cell track summaries
cell_summaries <- cells_corrected %>%
  group_by(cell_name) %>%
  summarise(
    n_steps        = n(),
    duration_min   = n_steps * FRAME_INTERVAL_MIN,
    mean_raw_speed = mean(raw_speed, na.rm = TRUE),
    mean_corr_speed = mean(corr_speed, na.rm = TRUE),
    median_corr_speed = median(corr_speed, na.rm = TRUE),
    total_raw_dist  = sum(sqrt(raw_dx^2 + raw_dy^2 + raw_dz^2), na.rm = TRUE),
    total_corr_dist = sum(sqrt(corr_dx^2 + corr_dy^2 + corr_dz^2), na.rm = TRUE),
    net_raw_disp    = sqrt(sum(raw_dx, na.rm = TRUE)^2 +
                           sum(raw_dy, na.rm = TRUE)^2 +
                           sum(raw_dz, na.rm = TRUE)^2),
    net_corr_disp   = sqrt(sum(corr_dx, na.rm = TRUE)^2 +
                           sum(corr_dy, na.rm = TRUE)^2 +
                           sum(corr_dz, na.rm = TRUE)^2),
    confinement_raw  = net_raw_disp / (total_raw_dist + 1e-10),
    confinement_corr = net_corr_disp / (total_corr_dist + 1e-10),
    .groups = "drop"
  )

cat(sprintf("\n  %d cells with corrected tracks\n", nrow(cell_summaries)))
cat(sprintf("  Median per-cell mean speed: %.3f µm/min (raw: %.3f)\n",
            median(cell_summaries$mean_corr_speed),
            median(cell_summaries$mean_raw_speed)))

# #############################################################################
# STEP 5: VISUALIZE CORRECTED MANUAL TRACKING
# #############################################################################

cat("\n=== STEP 5: CORRECTED TRACK VISUALIZATIONS ===\n")

# --- 5A: Raw vs corrected speed distributions ---
speed_compare <- bind_rows(
  cells_corrected %>% transmute(speed = raw_speed, type = "Raw (with drift)"),
  cells_corrected %>% transmute(speed = corr_speed, type = "Corrected (drift removed)")
)

p_speed_raw_corr <- ggplot(speed_compare, aes(x = speed, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Raw (with drift)" = "grey60",
                                "Corrected (drift removed)" = "#4DAF4A")) +
  coord_cartesian(xlim = c(0, quantile(speed_compare$speed, 0.99))) +
  labs(title = "Instantaneous Speed: Raw vs Drift-Corrected",
       x = "Speed (µm/min)", y = "Density", fill = NULL) +
  theme_manual()

# --- 5B: YSL verification — corrected YSL should be near zero ---
ysl_speed_compare <- bind_rows(
  ysl_corrected %>% transmute(speed = raw_speed,  type = "YSL raw"),
  ysl_corrected %>% transmute(speed = corr_speed, type = "YSL corrected")
)

p_ysl_verify <- ggplot(ysl_speed_compare, aes(x = speed, fill = type)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("YSL raw" = "#E66101", "YSL corrected" = "#4DAF4A")) +
  labs(title = "YSL Speed Verification",
       subtitle = "Corrected YSL ≈ 0 confirms drift removal works",
       x = "Speed (µm/min)", y = "Count", fill = NULL) +
  theme_manual()

# --- 5C: Corrected tracks — XY projection ---
# Reconstruct corrected positions by cumulating corrected displacements
corrected_positions <- cells_corrected %>%
  arrange(cell_name, FRAME) %>%
  group_by(cell_name) %>%
  mutate(
    corr_x = cumsum(corr_dx),
    corr_y = cumsum(corr_dy),
    corr_z = cumsum(corr_dz)
  ) %>%
  ungroup()

p_corr_xy <- ggplot(corrected_positions, aes(x = corr_x, y = corr_y,
                                              group = cell_name, color = cell_name)) +
  geom_path(linewidth = 0.5, alpha = 0.7) +
  geom_point(data = corrected_positions %>% group_by(cell_name) %>% slice(1),
             size = 2, shape = 1) +
  coord_fixed() +
  labs(title = "Corrected Cell Tracks (XY)",
       subtitle = "Embryo drift removed — true cell movement",
       x = "ΔX (µm)", y = "ΔY (µm)") +
  theme_manual() +
  theme(legend.position = "none")

# --- 5D: Speed over time ---
speed_over_time <- cells_corrected %>%
  mutate(time_bin = floor(time_min / 5) * 5) %>%
  group_by(time_bin) %>%
  summarise(
    mean_raw  = mean(raw_speed, na.rm = TRUE),
    mean_corr = mean(corr_speed, na.rm = TRUE),
    se_corr   = sd(corr_speed, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

p_speed_time <- ggplot() +
  geom_line(data = speed_over_time, aes(x = time_bin, y = mean_raw),
            color = "grey60", linewidth = 0.8, linetype = "dashed") +
  geom_ribbon(data = speed_over_time,
              aes(x = time_bin, ymin = mean_corr - se_corr,
                  ymax = mean_corr + se_corr),
              alpha = 0.3, fill = "#4DAF4A") +
  geom_line(data = speed_over_time, aes(x = time_bin, y = mean_corr),
            color = "#4DAF4A", linewidth = 0.8) +
  labs(title = "Mean Speed Over Time",
       subtitle = "Grey dashed = raw | Green = drift-corrected",
       x = "Time (min)", y = "Mean speed (µm/min)") +
  theme_manual()

# --- 5E: Per-cell mean speed bar plot ---
p_cell_speeds <- ggplot(cell_summaries %>%
                          arrange(mean_corr_speed) %>%
                          mutate(cell_name = factor(cell_name, levels = cell_name)),
                        aes(x = cell_name, y = mean_corr_speed)) +
  geom_col(fill = "#4DAF4A", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_corr_speed, ymax = mean_raw_speed),
                width = 0.3, color = "grey50", linewidth = 0.3) +
  coord_flip() +
  labs(title = "Per-Cell Mean Speed",
       subtitle = "Bar = corrected | Whisker extends to raw speed",
       x = NULL, y = "Mean speed (µm/min)") +
  theme_manual() +
  theme(axis.text.y = element_text(size = 6))

# --- 5F: Confinement ratio distribution ---
p_confinement <- ggplot(cell_summaries, aes(x = confinement_corr)) +
  geom_histogram(bins = 15, fill = "#4DAF4A", alpha = 0.7) +
  geom_vline(xintercept = c(0.3, 0.7), linetype = "dashed", color = "grey50") +
  labs(title = "Confinement Ratio (corrected)",
       subtitle = "<0.3 = confined; >0.7 = directed",
       x = "Confinement ratio", y = "Count") +
  theme_manual()

corr_page1 <- (p_speed_raw_corr + p_ysl_verify) / (p_corr_xy + p_speed_time) +
  plot_annotation(
    title = "Step 5: Drift-Corrected Manual Tracking",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "03_corrected_results.pdf"), corr_page1,
       width = 14, height = 12)
cat("  Saved: 03_corrected_results.pdf\n")

corr_page2 <- (p_cell_speeds + p_confinement) +
  plot_annotation(
    title = "Per-Cell Summaries (drift-corrected)",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "04_per_cell_summaries.pdf"), corr_page2,
       width = 14, height = 8)
cat("  Saved: 04_per_cell_summaries.pdf\n")

# --- Interactive 3D corrected tracks ---
fig_3d_manual <- plot_ly()
for (cn in unique(corrected_positions$cell_name)) {
  cdat <- corrected_positions %>% filter(cell_name == cn)
  fig_3d_manual <- fig_3d_manual %>%
    add_trace(data = cdat, x = ~corr_x, y = ~corr_y, z = ~corr_z,
              type = "scatter3d", mode = "lines",
              line = list(width = 3),
              name = cn, showlegend = TRUE,
              hoverinfo = "text",
              text = ~paste(cn, "| Frame:", FRAME))
}
fig_3d_manual <- fig_3d_manual %>%
  plotly::layout(
    title = list(text = "Manual Tracks — Drift-Corrected (3D)", font = list(size = 16)),
    scene = list(xaxis = list(title = "ΔX (µm)"),
                 yaxis = list(title = "ΔY (µm)"),
                 zaxis = list(title = "ΔZ (µm)"),
                 aspectmode = "data")
  )

htmlwidgets::saveWidget(fig_3d_manual,
                        file.path(OUTPUT_DIR, "05_interactive_3d_manual.html"),
                        selfcontained = TRUE)
cat("  Saved: 05_interactive_3d_manual.html\n")

# #############################################################################
# STEP 6: LOAD TRACKMATE & IMARIS SUMMARY STATS
# #############################################################################

cat("\n=== STEP 6: LOADING TRACKMATE & IMARIS DATA ===\n")

# --- TrackMate ---
cat("  Loading TrackMate...\n")
tm_spots_raw <- read_csv("4D_spots.csv", show_col_types = FALSE)[-c(1:3), ]
tm_tracks_raw <- read_csv("4D_tracks.csv", show_col_types = FALSE)[-c(1:3), ]

tm_spots <- tm_spots_raw %>%
  mutate(across(c(ID, TRACK_ID, POSITION_X, POSITION_Y, POSITION_Z, FRAME),
                ~as.numeric(.))) %>%
  filter(!is.na(TRACK_ID))

tm_tracks <- tm_tracks_raw %>%
  mutate(across(-LABEL, ~as.numeric(.))) %>%
  mutate(
    DURATION_MIN      = TRACK_DURATION * FRAME_INTERVAL_MIN,
    SPEED_MEAN_UM_MIN = TRACK_MEAN_SPEED * (60 / FRAME_INTERVAL_SEC),
    N_SPOTS           = NUMBER_SPOTS
  )

# QC (same as comparison script)
tm_tracks_clean <- tm_tracks %>%
  filter(N_SPOTS >= 5 & SPEED_MEAN_UM_MIN <= 10)
tm_spots_clean <- tm_spots %>%
  filter(TRACK_ID %in% tm_tracks_clean$TRACK_ID)

# Compute TM instantaneous speeds
tm_vel <- tm_spots_clean %>%
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

cat(sprintf("  TrackMate: %d tracks, %d spots after QC\n",
            nrow(tm_tracks_clean), nrow(tm_spots_clean)))

# --- Imaris ---
cat("  Loading Imaris...\n")
DATA_DIR_IM <- "all_tracks"
FILE_PREFIX <- "ome-tiff.companion-bin-8bit-crop"

read_imaris_csv <- function(stat_name) {
  path <- file.path(DATA_DIR_IM, paste0(FILE_PREFIX, "_", stat_name, ".csv"))
  read_csv(path, skip = 3, show_col_types = FALSE) %>%
    select(where(~ !all(is.na(.x))))
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

im_track_spd <- read_imaris_csv("Track_Speed_Mean") %>%
  rename(TRACK_SPEED_MEAN_RAW = `Track Speed Mean`, TRACK_ID = ID) %>%
  select(TRACK_ID, TRACK_SPEED_MEAN_RAW) %>%
  mutate(across(everything(), as.numeric))

im_track_str <- read_imaris_csv("Track_Straightness") %>%
  rename(CONFINEMENT_RATIO = `Track Straightness`, TRACK_ID = ID) %>%
  select(TRACK_ID, CONFINEMENT_RATIO) %>%
  mutate(across(everything(), as.numeric))

im_spots_per_track <- im_spots %>% count(TRACK_ID, name = "N_SPOTS")

im_tracks <- im_track_dur %>%
  left_join(im_track_spd, by = "TRACK_ID") %>%
  left_join(im_track_str, by = "TRACK_ID") %>%
  left_join(im_spots_per_track, by = "TRACK_ID") %>%
  mutate(
    DURATION_MIN      = TRACK_DURATION_FRAMES * FRAME_INTERVAL_SEC / 60,
    SPEED_MEAN_UM_MIN = TRACK_SPEED_MEAN_RAW / FRAME_INTERVAL_SEC * 60
  )

im_tracks_clean <- im_tracks %>%
  filter(N_SPOTS >= 5 & SPEED_MEAN_UM_MIN <= 10)
im_spots_clean <- im_spots %>%
  filter(TRACK_ID %in% im_tracks_clean$TRACK_ID)

# Compute Imaris instantaneous speeds
im_vel <- im_spots_clean %>%
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

cat(sprintf("  Imaris: %d tracks, %d spots after QC\n",
            nrow(im_tracks_clean), nrow(im_spots_clean)))

# #############################################################################
# STEP 7: 3-WAY COMPARISON
# #############################################################################

cat("\n=== STEP 7: 3-WAY COMPARISON ===\n")

# --- 7A: Instantaneous speed distributions ---
all_speeds_3way <- bind_rows(
  cells_corrected %>%
    transmute(inst_speed = corr_speed, software = "Manual (corrected)"),
  tm_vel %>%
    transmute(inst_speed = inst_speed, software = "TrackMate"),
  im_vel %>%
    transmute(inst_speed = inst_speed, software = "Imaris")
)

speed_summary_3way <- all_speeds_3way %>%
  group_by(software) %>%
  summarise(
    n           = n(),
    median_speed = median(inst_speed, na.rm = TRUE),
    mean_speed   = mean(inst_speed, na.rm = TRUE),
    sd_speed     = sd(inst_speed, na.rm = TRUE),
    q25          = quantile(inst_speed, 0.25, na.rm = TRUE),
    q75          = quantile(inst_speed, 0.75, na.rm = TRUE),
    q95          = quantile(inst_speed, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n  Instantaneous Speed Summary (µm/min):\n")
cat(sprintf("  %-22s  %8s  %8s  %8s  %8s  %8s\n",
            "Software", "Median", "Mean", "SD", "Q25", "Q75"))
for (i in 1:nrow(speed_summary_3way)) {
  cat(sprintf("  %-22s  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
              speed_summary_3way$software[i],
              speed_summary_3way$median_speed[i],
              speed_summary_3way$mean_speed[i],
              speed_summary_3way$sd_speed[i],
              speed_summary_3way$q25[i],
              speed_summary_3way$q75[i]))
}

# Density plot
p_3way_speed <- ggplot(all_speeds_3way, aes(x = inst_speed, fill = software)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = software_colors) +
  coord_cartesian(xlim = c(0, quantile(all_speeds_3way$inst_speed, 0.98, na.rm = TRUE))) +
  geom_vline(data = speed_summary_3way,
             aes(xintercept = median_speed, color = software),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  scale_color_manual(values = software_colors) +
  labs(title = "Instantaneous Speed — 3-Way Comparison",
       subtitle = "Dashed lines = medians",
       x = "Speed (µm/min)", y = "Density", fill = NULL) +
  theme_manual()

# Box plot
p_3way_box <- ggplot(all_speeds_3way, aes(x = software, y = inst_speed, fill = software)) +
  geom_boxplot(alpha = 0.6, outlier.alpha = 0.02, outlier.size = 0.3) +
  scale_fill_manual(values = software_colors) +
  coord_cartesian(ylim = c(0, quantile(all_speeds_3way$inst_speed, 0.98, na.rm = TRUE))) +
  labs(title = "Speed Distribution (box plot)",
       x = NULL, y = "Speed (µm/min)", fill = NULL) +
  theme_manual() +
  theme(legend.position = "none")

# --- 7B: Per-track mean speed comparison ---
tm_track_speeds <- tm_tracks_clean %>%
  transmute(mean_speed = SPEED_MEAN_UM_MIN, software = "TrackMate")
im_track_speeds <- im_tracks_clean %>%
  transmute(mean_speed = SPEED_MEAN_UM_MIN, software = "Imaris")
manual_track_speeds <- cell_summaries %>%
  transmute(mean_speed = mean_corr_speed, software = "Manual (corrected)")

all_track_speeds <- bind_rows(tm_track_speeds, im_track_speeds, manual_track_speeds)

p_3way_track_speed <- ggplot(all_track_speeds, aes(x = mean_speed, fill = software)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = software_colors) +
  coord_cartesian(xlim = c(0, quantile(all_track_speeds$mean_speed, 0.99, na.rm = TRUE))) +
  labs(title = "Per-Track Mean Speed — 3-Way",
       x = "Mean speed (µm/min)", y = "Density", fill = NULL) +
  theme_manual()

# --- 7C: Confinement ratio comparison ---
tm_confine <- tm_tracks_clean %>%
  transmute(confinement = CONFINEMENT_RATIO, software = "TrackMate")
im_confine <- im_tracks_clean %>%
  transmute(confinement = CONFINEMENT_RATIO, software = "Imaris")
manual_confine <- cell_summaries %>%
  transmute(confinement = confinement_corr, software = "Manual (corrected)")

all_confine <- bind_rows(tm_confine, im_confine, manual_confine)

p_3way_confine <- ggplot(all_confine, aes(x = confinement, fill = software)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = software_colors) +
  labs(title = "Confinement Ratio — 3-Way",
       x = "Confinement ratio", y = "Density", fill = NULL) +
  theme_manual()

# --- 7D: Track duration comparison ---
tm_dur <- tm_tracks_clean %>%
  transmute(duration_min = DURATION_MIN, software = "TrackMate")
im_dur <- im_tracks_clean %>%
  transmute(duration_min = DURATION_MIN, software = "Imaris")
manual_dur <- cell_summaries %>%
  transmute(duration_min = duration_min, software = "Manual (corrected)")

all_dur <- bind_rows(tm_dur, im_dur, manual_dur)

p_3way_dur <- ggplot(all_dur, aes(x = duration_min, fill = software)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = software_colors) +
  scale_x_log10() +
  labs(title = "Track Duration — 3-Way",
       x = "Duration (min, log)", y = "Density", fill = NULL) +
  theme_manual()

# Combine 3-way plots
threeway_combined <- (p_3way_speed + p_3way_box) / (p_3way_track_speed + p_3way_confine) / (p_3way_dur + plot_spacer()) +
  plot_annotation(
    title = "3-Way Comparison: Manual (ground truth) vs TrackMate vs Imaris",
    subtitle = sprintf("Manual: %d cells | TrackMate: %d tracks | Imaris: %d tracks",
                       nrow(cell_summaries), nrow(tm_tracks_clean), nrow(im_tracks_clean)),
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey50"))
  )

ggsave(file.path(OUTPUT_DIR, "06_3way_comparison.pdf"), threeway_combined,
       width = 14, height = 16)
cat("  Saved: 06_3way_comparison.pdf\n")

# --- 7E: Statistical tests ---
cat("\n  Pairwise Wilcoxon Tests (instantaneous speed):\n")

pairs <- list(
  c("Manual (corrected)", "TrackMate"),
  c("Manual (corrected)", "Imaris"),
  c("TrackMate", "Imaris")
)

for (pair in pairs) {
  v1 <- all_speeds_3way %>% filter(software == pair[1]) %>% pull(inst_speed)
  v2 <- all_speeds_3way %>% filter(software == pair[2]) %>% pull(inst_speed)
  wt <- wilcox.test(v1, v2, conf.int = FALSE)
  m1 <- median(v1, na.rm = TRUE)
  m2 <- median(v2, na.rm = TRUE)
  cat(sprintf("    %s vs %s: p = %.2e | medians: %.3f vs %.3f\n",
              pair[1], pair[2], wt$p.value, m1, m2))
}

# #############################################################################
# STEP 8: SUMMARY TABLE & EXPORT
# #############################################################################

cat("\n=== STEP 8: SUMMARY & EXPORT ===\n")

summary_table <- tibble(
  Metric = c("N tracks/cells",
             "Inst. speed median (µm/min)",
             "Inst. speed mean (µm/min)",
             "Per-track mean speed median",
             "Confinement ratio median",
             "Track duration median (min)"),
  Manual_corrected = c(
    nrow(cell_summaries),
    round(median(cells_corrected$corr_speed, na.rm = TRUE), 3),
    round(mean(cells_corrected$corr_speed, na.rm = TRUE), 3),
    round(median(cell_summaries$mean_corr_speed, na.rm = TRUE), 3),
    round(median(cell_summaries$confinement_corr, na.rm = TRUE), 3),
    round(median(cell_summaries$duration_min, na.rm = TRUE), 1)
  ),
  TrackMate = c(
    nrow(tm_tracks_clean),
    round(median(tm_vel$inst_speed, na.rm = TRUE), 3),
    round(mean(tm_vel$inst_speed, na.rm = TRUE), 3),
    round(median(tm_tracks_clean$SPEED_MEAN_UM_MIN, na.rm = TRUE), 3),
    round(median(tm_tracks_clean$CONFINEMENT_RATIO, na.rm = TRUE), 3),
    round(median(tm_tracks_clean$DURATION_MIN, na.rm = TRUE), 1)
  ),
  Imaris = c(
    nrow(im_tracks_clean),
    round(median(im_vel$inst_speed, na.rm = TRUE), 3),
    round(mean(im_vel$inst_speed, na.rm = TRUE), 3),
    round(median(im_tracks_clean$SPEED_MEAN_UM_MIN, na.rm = TRUE), 3),
    round(median(im_tracks_clean$CONFINEMENT_RATIO, na.rm = TRUE), 3),
    round(median(im_tracks_clean$DURATION_MIN, na.rm = TRUE), 1)
  )
)

cat("\n")
cat(strrep("=", 80), "\n")
cat("  3-WAY COMPARISON SUMMARY\n")
cat(strrep("=", 80), "\n\n")
print(summary_table, n = Inf, width = Inf)

# Export
write_csv(summary_table, file.path(OUTPUT_DIR, "3way_summary.csv"))
write_csv(cell_summaries, file.path(OUTPUT_DIR, "manual_cell_summaries.csv"))
write_csv(cells_corrected, file.path(OUTPUT_DIR, "manual_cells_corrected.csv"))
write_csv(drift_per_frame, file.path(OUTPUT_DIR, "ysl_drift_per_frame.csv"))
write_csv(speed_summary_3way, file.path(OUTPUT_DIR, "3way_speed_summary.csv"))

cat("\n  Exported CSV files to ", OUTPUT_DIR, "/\n")

# #############################################################################
# OUTPUT LIST
# #############################################################################

cat("\n")
cat(strrep("=", 80), "\n")
cat("  OUTPUT FILES (in ", OUTPUT_DIR, "/)\n")
cat(strrep("=", 80), "\n\n")
cat("  PDFs:\n")
cat("    01_raw_data_drift.pdf             — Raw tracks & YSL drift visualization\n")
cat("    02_drift_correction.pdf           — Drift magnitude, per-axis, cumulative path\n")
cat("    03_corrected_results.pdf          — Raw vs corrected speeds, YSL verification\n")
cat("    04_per_cell_summaries.pdf         — Per-cell speed bars, confinement histogram\n")
cat("    06_3way_comparison.pdf            — Speed, confinement, duration across 3 methods\n")
cat("\n  Interactive HTML:\n")
cat("    05_interactive_3d_manual.html     — 3D drift-corrected manual tracks\n")
cat("\n  Data:\n")
cat("    3way_summary.csv                  — Final summary table\n")
cat("    manual_cell_summaries.csv         — Per-cell summary statistics\n")
cat("    manual_cells_corrected.csv        — All corrected positions & velocities\n")
cat("    ysl_drift_per_frame.csv           — YSL drift correction vectors\n")
cat("    3way_speed_summary.csv            — Speed distribution summaries\n")

cat("\n")
cat(strrep("=", 80), "\n")
cat("  ANALYSIS COMPLETE\n")
cat(strrep("=", 80), "\n")
