# =============================================================================
# TrackMate 3D Cell Movement Analysis - Zebrafish/Medaka Embryo Nuclei Tracking
# =============================================================================
# Analysis pipeline for quality control, outlier detection, and characterization
# of nuclear movement during embryonic development
#
# BIOLOGICAL CONTEXT:
#   - Tracking nuclei in developing zebrafish/medaka embryos
#   - 3D imaging (SPIM) captures nuclear positions over time
#   - 1 frame = 30 seconds temporal resolution
#   - Nuclei move during cell division, migration, and morphogenesis
#   - Expected speeds: 0.5-5 µm/min for normal nuclear movement
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

FRAME_INTERVAL_SEC <- 30  # Each frame is 30 seconds
FRAME_INTERVAL_MIN <- 0.5 # 30 sec = 0.5 min

# =============================================================================
# COORDINATE TRANSFORMATION OPTIONS
# =============================================================================
#
# Microscopy coordinate systems may differ from image display conventions.
# Use these options to flip/mirror coordinates to match your original images.
#
# Set to TRUE to flip the corresponding axis:
#   - FLIP_X: mirrors left-right (multiplies X by -1)
#   - FLIP_Y: mirrors top-bottom (multiplies Y by -1)
#   - FLIP_Z: inverts depth (multiplies Z by -1)
#   - SWAP_XY: swaps X and Y axes (rotates 90°)
#
# Common scenarios:
#   - Image appears horizontally mirrored: set FLIP_X = TRUE
#   - Image appears vertically flipped: set FLIP_Y = TRUE
#   - Both mirrored and flipped: set FLIP_X = TRUE and FLIP_Y = TRUE
#

FLIP_X <- FALSE   # Set TRUE if image is horizontally mirrored
FLIP_Y <- FALSE   # Set TRUE if image is vertically flipped
FLIP_Z <- FALSE   # Set TRUE if Z-axis is inverted
SWAP_XY <- FALSE  # Set TRUE to swap X and Y axes


# Color palette for track categories
track_colors <- c(
  "High quality" = "#2166AC",
  "Medium quality" = "#B2182B",
  "Low quality" = "#762A83",
  "Filtered" = "grey85"
)

# Movement type colors
movement_colors <- c(
  "Directed" = "#2166AC",
  "Confined" = "#B2182B",
  "Random" = "#762A83",
  "NS" = "grey85"
)

# Cluster colors (using viridis for many categories)
cluster_pal <- viridis::viridis(6)

# Base plotting theme
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
# PART 1: DATA LOADING AND PREPROCESSING
# =============================================================================
#
# INTERPRETATION:
#   TrackMate outputs spots (individual detections) and tracks (linked spots).
#   Each spot is a nucleus detected at a specific (X, Y, Z, T) position.
#   Tracks connect the same nucleus across time frames.
#
#   Quality control is essential because:
#   - False detections (noise, debris) create spurious short tracks
#   - Tracking errors (ID swaps) create unrealistically fast movements
#   - Dividing cells may cause track splits/merges
#

cat("\n=== PART 1: DATA LOADING ===\n")

# Load data - skipping TrackMate metadata rows
spots_raw <- read_csv("oriented_spots.csv", show_col_types = FALSE)[-c(1:3),]
tracks_raw <- read_csv("oriented_tracks.csv", show_col_types = FALSE)[-c(1:3),]

# Convert to numeric
spots <- spots_raw %>%
  mutate(across(c(ID, TRACK_ID, QUALITY, POSITION_X, POSITION_Y, POSITION_Z,
                  POSITION_T, FRAME, RADIUS, MEAN_INTENSITY_CH1,
                  MEDIAN_INTENSITY_CH1, MIN_INTENSITY_CH1, MAX_INTENSITY_CH1,
                  TOTAL_INTENSITY_CH1, STD_INTENSITY_CH1, CONTRAST_CH1, SNR_CH1),
                ~as.numeric(.))) %>%
  filter(!is.na(TRACK_ID))

tracks <- tracks_raw %>%
  mutate(across(-LABEL, ~as.numeric(.)))

# =============================================================================
# Apply coordinate transformations if requested
# =============================================================================

if (FLIP_X || FLIP_Y || FLIP_Z || SWAP_XY) {
  cat("Applying coordinate transformations:\n")
  if (FLIP_X) cat("  - Flipping X axis (horizontal mirror)\n")
  if (FLIP_Y) cat("  - Flipping Y axis (vertical flip)\n")
  if (FLIP_Z) cat("  - Flipping Z axis (depth inversion)\n")
  if (SWAP_XY) cat("  - Swapping X and Y axes\n")

  # Store original coordinates for reference
  spots <- spots %>%
    mutate(
      POSITION_X_ORIG = POSITION_X,
      POSITION_Y_ORIG = POSITION_Y,
      POSITION_Z_ORIG = POSITION_Z
    )

  # Apply transformations
  if (SWAP_XY) {
    spots <- spots %>%
      mutate(
        temp_x = POSITION_X,
        POSITION_X = POSITION_Y,
        POSITION_Y = temp_x
      ) %>%
      select(-temp_x)
  }

  if (FLIP_X) spots <- spots %>% mutate(POSITION_X = -POSITION_X)
  if (FLIP_Y) spots <- spots %>% mutate(POSITION_Y = -POSITION_Y)
  if (FLIP_Z) spots <- spots %>% mutate(POSITION_Z = -POSITION_Z)

  # Shift to positive coordinates (optional, for cleaner visualization)
  spots <- spots %>%
    mutate(
      POSITION_X = POSITION_X - min(POSITION_X, na.rm = TRUE),
      POSITION_Y = POSITION_Y - min(POSITION_Y, na.rm = TRUE),
      POSITION_Z = POSITION_Z - min(POSITION_Z, na.rm = TRUE)
    )

  cat("\n")
}

# =============================================================================
# CONVERT TIME UNITS: frames → real time
# =============================================================================

tracks <- tracks %>%
  mutate(
    # Convert duration from frames to minutes
    DURATION_MIN = TRACK_DURATION * FRAME_INTERVAL_MIN,

    # Convert speeds from µm/frame to µm/min
    # Original: distance moved per frame; multiply by frames/min
    SPEED_MEAN_UM_MIN = TRACK_MEAN_SPEED * (60 / FRAME_INTERVAL_SEC),
    SPEED_MAX_UM_MIN = TRACK_MAX_SPEED * (60 / FRAME_INTERVAL_SEC),
    SPEED_MEDIAN_UM_MIN = TRACK_MEDIAN_SPEED * (60 / FRAME_INTERVAL_SEC)
  )

cat("Data Summary:\n")
cat("  Total tracks:", nrow(tracks), "\n")
cat("  Total spots:", nrow(spots), "\n")
cat("  Time covered:", max(spots$FRAME, na.rm = TRUE) * FRAME_INTERVAL_MIN, "min\n")
cat("  Frame interval:", FRAME_INTERVAL_SEC, "sec\n\n")

# =============================================================================
# PART 2: DATA EXPLORATION
# =============================================================================
#
# INTERPRETATION:
#   Understanding the distribution of track properties helps identify:
#   - Normal range of nuclear movement speeds
#   - Typical track durations (how long nuclei can be followed)
#   - Potential tracking artifacts (outliers)
#

cat("=== PART 2: DATA EXPLORATION ===\n\n")

# Summary with converted units
track_summary <- tracks %>%
  summarise(
    # Duration in minutes
    duration_min_median = median(DURATION_MIN, na.rm = TRUE),
    duration_min_mean = mean(DURATION_MIN, na.rm = TRUE),
    duration_min_max = max(DURATION_MIN, na.rm = TRUE),

    # Speed in µm/min
    speed_um_min_median = median(SPEED_MEAN_UM_MIN, na.rm = TRUE),
    speed_um_min_mean = mean(SPEED_MEAN_UM_MIN, na.rm = TRUE),
    speed_um_min_max = max(SPEED_MAX_UM_MIN, na.rm = TRUE),

    # Displacement and distance
    displacement_median = median(TRACK_DISPLACEMENT, na.rm = TRUE),
    total_dist_median = median(TOTAL_DISTANCE_TRAVELED, na.rm = TRUE),

    # Confinement
    confinement_median = median(CONFINEMENT_RATIO, na.rm = TRUE),

    # Track length
    n_spots_median = median(NUMBER_SPOTS, na.rm = TRUE)
  )

cat("Track Statistics (converted to real time):\n")
cat(sprintf("  Median duration: %.1f min (%.0f frames)\n",
            track_summary$duration_min_median,
            track_summary$duration_min_median / FRAME_INTERVAL_MIN))
cat(sprintf("  Median speed: %.2f µm/min\n", track_summary$speed_um_min_median))
cat(sprintf("  Median displacement: %.1f µm\n", track_summary$displacement_median))
cat(sprintf("  Median confinement ratio: %.2f\n", track_summary$confinement_median))
cat(sprintf("  Median track length: %.0f spots\n\n", track_summary$n_spots_median))

# =============================================================================
# PART 3: QUALITY CONTROL AND FILTERING
# =============================================================================
#
# BIOLOGICAL RATIONALE FOR FILTERS:
#
# 1. MINIMUM TRACK LENGTH (≥5 spots = ≥2.5 min):
#    - Very short tracks are likely noise or transient detections
#    - Need several timepoints to calculate meaningful velocities
#    - 2.5 min is enough to capture directed movement vs random diffusion
#
# 2. MAXIMUM SPEED (≤10 µm/min):
#    - Even fast-migrating cells rarely exceed 10 µm/min
#    - Higher speeds usually indicate tracking errors (ID swaps)
#    - Nuclear "jumps" between frames = linking errors
#
# 3. MINIMUM DURATION (≥1.5 min = 3 frames):
#    - Single-frame "tracks" have no velocity information
#    - Need multiple timepoints to assess movement pattern
#
# 4. QUALITY SCORE (TrackMate's detection confidence):
#    - Low quality = dim/blurry nuclei, less reliable detection
#    - Keep top 90% by quality
#

cat("=== PART 3: QUALITY CONTROL ===\n\n")




# Define biologically sensible thresholds
MIN_SPOTS <- 5
MIN_DURATION_MIN <- 3
QUALITY_PERCENTILE <- 0.10  # Remove bottom 10% quality
MEAN_SPEED_THRESHOLD <- 5    # µm/min - tracks with MEAN above this are likely errors
MEDIAN_SPEED_THRESHOLD <- 5   # µm/min - more robust estimate
FAST_FRACTION_THRESHOLD <- 0.2 # Remove if >50% of steps are too fast
STEP_SPEED_THRESHOLD <- 5     # µm/min - what counts as a "fast" step

quality_threshold <- quantile(tracks$TRACK_MEAN_QUALITY, QUALITY_PERCENTILE, na.rm = TRUE)

# First, calculate per-track speed statistics from spots
# This gives us finer control than just using TRACK_MAX_SPEED
spot_speed_stats <- spots %>%
  arrange(TRACK_ID, FRAME) %>%
  group_by(TRACK_ID) %>%
  mutate(
    dx = POSITION_X - lag(POSITION_X),
    dy = POSITION_Y - lag(POSITION_Y),
    dz = POSITION_Z - lag(POSITION_Z),
    dt_frames = FRAME - lag(FRAME),
    step_dist = sqrt(dx^2 + dy^2 + dz^2),
    step_speed_um_min = (step_dist / dt_frames) * (60 / FRAME_INTERVAL_SEC)
  ) %>%
  filter(!is.na(step_speed_um_min) & dt_frames > 0) %>%
  summarise(
    n_steps = n(),
    mean_step_speed = mean(step_speed_um_min, na.rm = TRUE),
    median_step_speed = median(step_speed_um_min, na.rm = TRUE),
    max_step_speed = max(step_speed_um_min, na.rm = TRUE),
    sd_step_speed = sd(step_speed_um_min, na.rm = TRUE),
    pct_fast_steps = mean(step_speed_um_min > STEP_SPEED_THRESHOLD, na.rm = TRUE),
    .groups = "drop"
  )

# Merge speed stats back to tracks
tracks <- tracks %>%
  left_join(spot_speed_stats, by = "TRACK_ID")

# Flag problematic tracks
tracks <- tracks %>%
  mutate(
    # Individual flags
    flag_too_short = NUMBER_SPOTS < MIN_SPOTS,
    flag_too_brief = DURATION_MIN < MIN_DURATION_MIN,

    # NEW: Consistently too fast (not just one fast moment)
    flag_consistently_fast = (
      mean_step_speed > MEAN_SPEED_THRESHOLD |
      median_step_speed > MEDIAN_SPEED_THRESHOLD |
      pct_fast_steps > FAST_FRACTION_THRESHOLD
    ),
    # Replace NA with FALSE for tracks with too few steps
    flag_consistently_fast = ifelse(is.na(flag_consistently_fast), FALSE, flag_consistently_fast),

    flag_low_quality = TRACK_MEAN_QUALITY < quality_threshold,

    # Summary
    n_flags = flag_too_short + flag_too_brief + flag_consistently_fast + flag_low_quality,
    passes_qc = n_flags == 0
  )

# Summary of flags
cat("Quality Control Flags:\n")
cat(sprintf("  Too short (<%d spots): %d tracks (%.1f%%)\n",
            MIN_SPOTS, sum(tracks$flag_too_short), 100 * mean(tracks$flag_too_short)))
cat(sprintf("  Too brief (<%.1f min): %d tracks (%.1f%%)\n",
            MIN_DURATION_MIN, sum(tracks$flag_too_brief), 100 * mean(tracks$flag_too_brief)))
cat(sprintf("  Consistently too fast: %d tracks (%.1f%%)\n",
            sum(tracks$flag_consistently_fast, na.rm = TRUE),
            100 * mean(tracks$flag_consistently_fast, na.rm = TRUE)))
cat("    (mean >%.0f OR median >%.0f OR >%.0f%% steps >%.0f µm/min)\n",
    MEAN_SPEED_THRESHOLD, MEDIAN_SPEED_THRESHOLD,
    100*FAST_FRACTION_THRESHOLD, STEP_SPEED_THRESHOLD)
cat(sprintf("  Low quality (bottom 10%%): %d tracks (%.1f%%)\n",
            sum(tracks$flag_low_quality), 100 * mean(tracks$flag_low_quality)))

# Show what we're keeping vs filtering for fast tracks specifically
cat("\nSpeed filtering details:\n")
cat(sprintf("  Tracks with occasional fast burst (kept): %d\n",
            sum(!tracks$flag_consistently_fast & tracks$max_step_speed > STEP_SPEED_THRESHOLD, na.rm = TRUE)))
cat(sprintf("  Tracks consistently fast (removed): %d\n",
            sum(tracks$flag_consistently_fast, na.rm = TRUE)))

# Apply filters
tracks_clean <- tracks %>% filter(passes_qc)
spots_clean <- spots %>% filter(TRACK_ID %in% tracks_clean$TRACK_ID)

cat(sprintf("\nFiltering Results:\n"))
cat(sprintf("  Original: %d tracks, %d spots\n", nrow(tracks), nrow(spots)))
cat(sprintf("  Filtered: %d tracks (%.1f%%), %d spots\n",
            nrow(tracks_clean), 100 * nrow(tracks_clean) / nrow(tracks), nrow(spots_clean)))

# =============================================================================
# PART 4: QUALITY CONTROL VISUALIZATIONS
# =============================================================================
#
# INTERPRETATION:
#   These plots help verify that our filtering is sensible:
#   - Distribution shifts show what we're removing
#   - Outliers in speed plot = tracking errors
#   - Duration vs spots should be linear (lost tracks = deviations)
#

cat("\n=== PART 4: QC VISUALIZATIONS ===\n")

if (!dir.exists("analysis_output")) dir.create("analysis_output")

# Plot 1: Track length distribution (before/after)
p1 <- ggplot() +
  geom_histogram(data = tracks, aes(x = NUMBER_SPOTS, fill = "Before filtering"),
                 alpha = 0.5, bins = 50) +
  geom_histogram(data = tracks_clean, aes(x = NUMBER_SPOTS, fill = "After filtering"),
                 alpha = 0.5, bins = 50) +
  geom_vline(xintercept = MIN_SPOTS, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_fill_manual(values = c("Before filtering" = "grey50", "After filtering" = "#2166AC")) +
  scale_x_log10() +
  labs(
    title = "Track Length Distribution",
    subtitle = sprintf("Dashed line: minimum %d spots (%.1f min)", MIN_SPOTS, MIN_SPOTS * FRAME_INTERVAL_MIN),
    x = "Number of spots per track (log scale)",
    y = "Count",
    fill = NULL
  ) +
  theme_track() +
  theme(legend.position = "top")

# Plot 2: Speed distribution
#
# INTERPRETATION:
#   - Main peak = typical nuclear movement speed
#   - Long tail = tracking errors or fast migrations
#   - Compare to known nuclear speeds (0.5-5 µm/min typical)
#
p2 <- ggplot(tracks_clean, aes(x = SPEED_MEAN_UM_MIN)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "#2166AC", alpha = 0.7) +
  geom_density(color = "#B2182B", linewidth = 1) +
  geom_vline(xintercept = median(tracks_clean$SPEED_MEAN_UM_MIN, na.rm = TRUE),
             linetype = "dashed", color = "grey30", linewidth = 0.4) +
  labs(
    title = "Nuclear Migration Speed Distribution (3D)",
    subtitle = sprintf("Median: %.2f µm/min (XYZ)",
                      median(tracks_clean$SPEED_MEAN_UM_MIN, na.rm = TRUE)),
    x = "Mean speed (µm/min)",
    y = "Density"
  ) +
  theme_track()

# Plot 3: Duration vs Max Speed (identify outliers)
#
# INTERPRETATION:
#   - Most tracks cluster at moderate speeds
#   - Points in upper region = suspiciously fast (tracking errors)
#   - Filtering removes these outliers
#
p3 <- ggplot(tracks, aes(x = DURATION_MIN, y = SPEED_MAX_UM_MIN)) +
  geom_point(aes(color = ifelse(passes_qc, "Kept", "Filtered")),
             alpha = 0.3, size = 0.5) +
  geom_hline(yintercept = MEAN_SPEED_THRESHOLD, linetype = "dashed", color = "#B2182B", linewidth = 0.3) +
  geom_vline(xintercept = MIN_DURATION_MIN, linetype = "dashed", color = "#B2182B", linewidth = 0.3) +
  scale_color_manual(values = c("Kept" = "#2166AC", "Filtered" = "grey70"),
                     breaks = c("Kept", "Filtered")) +
  scale_x_log10() +
  labs(
    title = "Track Duration vs Maximum Speed",
    subtitle = "Red dashed lines: filtering thresholds",
    x = "Duration (min, log scale)",
    y = "Max speed (µm/min)",
    color = NULL
  ) +
  theme_track()

# Plot 4: Quality score distribution
p4 <- ggplot(tracks, aes(x = TRACK_MEAN_QUALITY, fill = passes_qc)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = quality_threshold, linetype = "dashed", color = "#B2182B", linewidth = 0.3) +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "grey70"),
                    labels = c("TRUE" = "Kept", "FALSE" = "Filtered")) +
  labs(
    title = "Detection Quality Score",
    subtitle = "TrackMate confidence in nucleus detection",
    x = "Mean quality score",
    y = "Count",
    fill = NULL
  ) +
  theme_track()

# Combine QC plots
qc_combined <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Quality Control Summary",
    subtitle = sprintf("Kept %d of %d tracks (%.1f%%)",
                      nrow(tracks_clean), nrow(tracks), 100 * nrow(tracks_clean) / nrow(tracks)),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
    )
  )

ggsave("analysis_output/01_quality_control.pdf", qc_combined, width = 10, height = 10)
cat("  Saved: 01_quality_control.pdf\n")

# =============================================================================
# PART 5: VELOCITY ANALYSIS
# =============================================================================
#
# BIOLOGICAL INTERPRETATION:
#   Nuclear velocities reflect cell behavior:
#   - Constant speed → directed migration (e.g., during gastrulation)
#   - Variable speed → random movement or confined exploration
#   - Speed changes over time → developmental transitions
#
#   We calculate instantaneous velocities from spot-to-spot displacements.
#

cat("\n=== PART 5: VELOCITY ANALYSIS ===\n")

# Calculate instantaneous velocities
spots_velocity <- spots_clean %>%
  arrange(TRACK_ID, FRAME) %>%
  group_by(TRACK_ID) %>%
  mutate(
    # Time difference (in minutes)
    dt_min = (FRAME - lag(FRAME)) * FRAME_INTERVAL_MIN,

    # Displacements
    dx = POSITION_X - lag(POSITION_X),
    dy = POSITION_Y - lag(POSITION_Y),
    dz = POSITION_Z - lag(POSITION_Z),
    displacement = sqrt(dx^2 + dy^2 + dz^2),

    # Instantaneous speed (µm/min)
    inst_speed_um_min = displacement / dt_min,

    # Velocity components (µm/min)
    vx = dx / dt_min,
    vy = dy / dt_min,
    vz = dz / dt_min,

    # Time from track start (minutes)
    time_min = (FRAME - min(FRAME)) * FRAME_INTERVAL_MIN
  ) %>%
  ungroup()

# Speed statistics
speed_stats <- spots_velocity %>%
  filter(!is.na(inst_speed_um_min) ) %>%  # Remove outliers
  summarise(
    mean = mean(inst_speed_um_min),
    median = median(inst_speed_um_min),
    sd = sd(inst_speed_um_min),
    q95 = quantile(inst_speed_um_min, 0.95)
  )

cat(sprintf("Instantaneous Speed Statistics:\n"))
cat(sprintf("  Mean: %.2f µm/min\n", speed_stats$mean))
cat(sprintf("  Median: %.2f µm/min\n", speed_stats$median))
cat(sprintf("  95th percentile: %.2f µm/min\n", speed_stats$q95))

# Plot 5: Instantaneous speed distribution
#
# INTERPRETATION:
#   - Peak location = typical nuclear speed
#   - Width = variability in movement
#   - Tail = fast movements (division, migration)
#
p5 <- ggplot(spots_velocity %>% filter(!is.na(inst_speed_um_min)),
             aes(x = inst_speed_um_min)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "#2166AC", alpha = 0.7) +
  geom_density(color = "#B2182B", linewidth = 1) +
  geom_vline(xintercept = speed_stats$median, linetype = "dashed", color = "grey30", linewidth = 0.4) +
  labs(
    title = "Instantaneous Nuclear Speed (3D)",
    subtitle = sprintf("Median: %.2f µm/min | 95th pct: %.2f µm/min (XYZ displacement)",
                       speed_stats$median, speed_stats$q95),
    x = "Speed (µm/min)",
    y = "Density"
  ) +
  theme_track()

# Plot 6: Speed over developmental time (10-min bins)
speed_over_time <- spots_velocity %>%
  filter(!is.na(inst_speed_um_min)) %>%
  mutate(time_bin_min = floor(FRAME * FRAME_INTERVAL_MIN / 10) * 10) %>%
  group_by(time_bin_min) %>%
  summarise(
    mean_speed = mean(inst_speed_um_min, na.rm = TRUE),
    se_speed = sd(inst_speed_um_min, na.rm = TRUE) / sqrt(n()),
    n = n()
  )

p6 <- ggplot(speed_over_time, aes(x = time_bin_min, y = mean_speed)) +
  geom_ribbon(aes(ymin = mean_speed - se_speed, ymax = mean_speed + se_speed),
              alpha = 0.3, fill = "#2166AC") +
  geom_line(color = "#2166AC", linewidth = 1) +
  labs(
    title = "Mean Instantenous Speed Over Time",
    subtitle = "Shaded: SE | Speed computed from XYZ displacements",
    x = "Time (min from start)",
    y = "Mean speed (µm/min)"
  ) +
  theme_track()

# Combine velocity plots
velocity_combined <- (p5 + p6) +
  plot_annotation(
    title = "Nuclear Velocity Analysis",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

ggsave("analysis_output/02_velocity_analysis.pdf", velocity_combined, width = 18, height = 5)
cat("  Saved: 02_velocity_analysis.pdf\n")

# =============================================================================
# PART 6: DIRECTIONALITY ANALYSIS
# =============================================================================
#
# BIOLOGICAL INTERPRETATION:
#   Nuclear movement directionality reveals:
#   - Random movement (high turning angles) → confined/exploratory behavior
#   - Persistent movement (low turning angles) → directed migration
#   - Preferred directions → tissue-level organization forces
#
#   Turning angle: angle between consecutive displacement vectors
#   - 0° = straight line (persistent)
#   - 90° = random (no correlation)
#   - 180° = reversal (oscillatory)
#

cat("\n=== PART 6: DIRECTIONALITY ANALYSIS ===\n")

# Calculate turning angles
spots_direction <- spots_velocity %>%
  filter(!is.na(dx) & !is.na(dy) & !is.na(dz)) %>%
  group_by(TRACK_ID) %>%
  mutate(
    # Previous displacement vector
    prev_dx = lag(dx),
    prev_dy = lag(dy),
    prev_dz = lag(dz),

    # Dot product and magnitudes
    dot_product = dx * prev_dx + dy * prev_dy + dz * prev_dz,
    mag_current = sqrt(dx^2 + dy^2 + dz^2),
    mag_prev = sqrt(prev_dx^2 + prev_dy^2 + prev_dz^2),

    # Turning angle (radians → degrees)
    turning_angle = acos(pmin(1, pmax(-1, dot_product / (mag_current * mag_prev + 1e-10)))) * 180 / pi,

    # XY direction (azimuthal angle)
    direction_xy = atan2(dy, dx) * 180 / pi
  ) %>%
  ungroup()

# Turning angle summary
turning_stats <- spots_direction %>%
  filter(!is.na(turning_angle)) %>%
  summarise(
    mean_angle = mean(turning_angle),
    median_angle = median(turning_angle)
  )

cat(sprintf("Turning Angle Statistics:\n"))
cat(sprintf("  Mean: %.1f°\n", turning_stats$mean_angle))
cat(sprintf("  Median: %.1f°\n", turning_stats$median_angle))
cat("  (90° = random; <90° = persistent; >90° = reversing)\n")

# Plot 7: Turning angle distribution
#
# INTERPRETATION:
#   - Peak near 90° → random/diffusive movement
#   - Peak below 90° → directed/persistent movement
#   - Bimodal → mixed populations (some directed, some random)
#
p7 <- ggplot(spots_direction %>% filter(!is.na(turning_angle)),
             aes(x = turning_angle)) +
  geom_histogram(aes(y = after_stat(density)), bins = 36, fill = "#2166AC", alpha = 0.7) +
  geom_vline(xintercept = 90, linetype = "dashed", color = "#B2182B", linewidth = 0.4) +
  labs(
    title = "Distribution of Turning Angles",
    subtitle = "90° = random; <90° = persistent; >90° = reversing",
    x = "Turning angle (degrees)",
    y = "Density"
  ) +
  theme_track()

# Plot 8: Rose diagram of XY directions
#
# INTERPRETATION:
#   - Uniform → no preferred direction (isotropic movement)
#   - Peaks → preferred migration axes (tissue polarity)
#
p8 <- ggplot(spots_direction %>% filter(!is.na(direction_xy)),
             aes(x = direction_xy)) +
  geom_histogram(bins = 36, fill = "#2166AC", color = "white", alpha = 0.7) +
  coord_polar(start = 0) +
  scale_x_continuous(limits = c(-180, 180), breaks = seq(-180, 135, 45)) +
  labs(
    title = "Movement Direction (XY plane)",
    subtitle = "Rose diagram of migration angles",
    x = "Direction (degrees)",
    y = "Count"
  ) +
  theme_track() +
  theme(axis.text.y = element_blank())

# Combine direction plots
direction_combined <- (p7 + p8) +
  plot_annotation(
    title = "Directionality Analysis",
    subtitle = "How nuclei change direction during movement",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
    )
  )

ggsave("analysis_output/03_directionality.pdf", direction_combined, width = 12, height = 5)
cat("  Saved: 03_directionality.pdf\n")

# =============================================================================
# PART 7: SPATIAL VELOCITY PATTERNS
# =============================================================================
#
# BIOLOGICAL INTERPRETATION:
#   Spatial patterns in velocity reveal:
#   - Morphogenetic flows (collective cell movements)
#   - Regional differences (e.g., faster at embryo margins)
#   - Tissue organization (convergent/divergent flows)
#

cat("\n=== PART 7: SPATIAL VELOCITY PATTERNS ===\n")

# Bin space for local velocity statistics
spatial_velocity <- spots_velocity %>%
  filter(!is.na(inst_speed_um_min)) %>%
  mutate(
    x_bin = round(POSITION_X / 15) * 15,  # 15 µm bins
    y_bin = round(POSITION_Y / 15) * 15,
    z_bin = round(POSITION_Z / 15) * 15
  )

# Local velocity statistics (XY projection)
local_stats_xy <- spatial_velocity %>%
  group_by(x_bin, y_bin) %>%
  summarise(
    mean_speed = mean(inst_speed_um_min, na.rm = TRUE),
    mean_vx = mean(vx, na.rm = TRUE),
    mean_vy = mean(vy, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 20)

# Plot 9: Spatial speed heatmap
#
# INTERPRETATION:
#   - Hot spots = regions of high nuclear motility
#   - Cold regions = quiescent or stationary nuclei
#   - Gradients = transitional zones
#
p9 <- ggplot(local_stats_xy, aes(x = x_bin, y = y_bin, fill = mean_speed)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Speed\n(µm/min)", option = "plasma") +
  coord_fixed() +
  labs(
    title = "Spatial Distribution of Nuclear Speed",
    subtitle = "XY projection; color = mean local speed",
    x = "X position (µm)",
    y = "Y position (µm)"
  ) +
  theme_track()

# Plot 10: Velocity vector field
#
# INTERPRETATION:
#   - Arrows show local flow direction
#   - Aligned arrows = collective migration
#   - Random arrows = disorganized movement
#
p10 <- ggplot(local_stats_xy %>% filter(n >= 50),
              aes(x = x_bin, y = y_bin)) +
  geom_segment(aes(xend = x_bin + mean_vx * 5, yend = y_bin + mean_vy * 5,
                   color = mean_speed),
               arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.4) +
  scale_color_viridis_c(name = "Speed\n(µm/min)", option = "plasma") +
  coord_fixed() +
  labs(
    title = "Local Velocity Vector Field",
    subtitle = "Arrows show mean flow direction (scaled 5x)",
    x = "X position (µm)",
    y = "Y position (µm)"
  ) +
  theme_track()

# Combine spatial plots
spatial_combined <- (p9 + p10) +
  plot_annotation(
    title = "Spatial Velocity Patterns (15 um clusters)",
    subtitle = "Where nuclei move faster and in what direction",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
    )
  )

ggsave("analysis_output/04_spatial_velocity.pdf", spatial_combined, width = 12, height = 6)
cat("  Saved: 04_spatial_velocity.pdf\n")

# =============================================================================
# PART 8: MEAN SQUARED DISPLACEMENT (MSD) ANALYSIS
# =============================================================================
#
# BIOLOGICAL INTERPRETATION:
#   MSD reveals the TYPE of movement:
#   - MSD ~ t^α where α (alpha) indicates movement mode:
#     α ≈ 1: Normal diffusion (random walk)
#     α < 1: Subdiffusion (confined, e.g., in crowded tissue)
#     α > 1: Superdiffusion (directed migration)
#
#   This is key for understanding whether nuclei are:
#   - Freely diffusing (random)
#   - Confined by surrounding cells
#   - Actively migrating with a bias
#

cat("\n=== PART 8: MSD ANALYSIS ===\n")

# Calculate MSD for each track
calculate_msd <- function(df, max_lag = 30) {
  n <- nrow(df)
  if (n < MIN_SPOTS) return(NULL)

  max_lag <- min(max_lag, n - 1)

  msd_values <- sapply(1:max_lag, function(lag) {
    displacements <- (df$POSITION_X[(lag+1):n] - df$POSITION_X[1:(n-lag)])^2 +
      (df$POSITION_Y[(lag+1):n] - df$POSITION_Y[1:(n-lag)])^2 +
      (df$POSITION_Z[(lag+1):n] - df$POSITION_Z[1:(n-lag)])^2
    mean(displacements, na.rm = TRUE)
  })

  tibble(
    lag_frames = 1:max_lag,
    lag_min = 1:max_lag * FRAME_INTERVAL_MIN,
    msd = msd_values
  )
}

# Sample tracks for MSD (memory efficient)
set.seed(42)
sample_tracks <- sample(unique(spots_clean$TRACK_ID), min(10000, n_distinct(spots_clean$TRACK_ID)))

msd_data <- spots_clean %>%
  filter(TRACK_ID %in% sample_tracks) %>%
  arrange(TRACK_ID, FRAME) %>%
  group_by(TRACK_ID) %>%
  group_modify(~{
    result <- calculate_msd(.x)
    if (is.null(result)) return(tibble())
    result
  }) %>%
  ungroup()

# Ensemble MSD
ensemble_msd <- msd_data %>%
  group_by(lag_min) %>%
  summarise(
    mean_msd = mean(msd, na.rm = TRUE),
    se_msd = sd(msd, na.rm = TRUE) / sqrt(n()),
    n = n()
  ) %>%
  filter(n >= 50)

# Fit power law: MSD = D * t^alpha
msd_fit <- lm(log(mean_msd) ~ log(lag_min), data = ensemble_msd %>% filter(lag_min <= 10))
alpha <- coef(msd_fit)[2]
D_apparent <- exp(coef(msd_fit)[1]) / 6  # 3D diffusion coefficient

cat(sprintf("MSD Analysis Results:\n"))
cat(sprintf("  Power law exponent (α): %.3f\n", alpha))
cat(sprintf("  Interpretation: %s\n",
            ifelse(alpha < 0.8, "Subdiffusive (confined movement)",
                   ifelse(alpha > 1.2, "Superdiffusive (directed migration)",
                          "Normal diffusion (random walk)"))))
cat(sprintf("  Apparent diffusion coefficient: %.3f µm²/min\n", D_apparent))

# Plot 11: MSD curves
#
# INTERPRETATION:
#   - Individual gray lines = single track MSDs (variability)
#   - Red line = ensemble average (population behavior)
#   - Slope on log-log plot = α (movement type)
#
p11 <- ggplot() +
  geom_line(data = msd_data %>% filter(TRACK_ID %in% sample(unique(TRACK_ID), 50)),
            aes(x = lag_min, y = msd, group = TRACK_ID), alpha = 0.2, color = "grey50") +
  geom_ribbon(data = ensemble_msd,
              aes(x = lag_min, ymin = mean_msd - se_msd, ymax = mean_msd + se_msd),
              alpha = 0.3, fill = "#B2182B") +
  geom_line(data = ensemble_msd, aes(x = lag_min, y = mean_msd),
            color = "#B2182B", linewidth = 1.2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Mean Squared Displacement Analysis",
    subtitle = sprintf("alpha = %.2f (%s)", alpha,
                      ifelse(alpha < 0.8, "subdiffusive",
                             ifelse(alpha > 1.2, "superdiffusive", "diffusive"))),
    x = "Time lag (min, log scale)",
    y = "MSD (µm², log scale)"
  ) +
  theme_track()

ggsave("analysis_output/05_msd_analysis.pdf", p11, width = 8, height = 6)
cat("  Saved: 05_msd_analysis.pdf\n")

# =============================================================================
# PART 9: TRACK CLASSIFICATION BY MOVEMENT TYPE
# =============================================================================
#
# BIOLOGICAL INTERPRETATION:
#   Different nuclei show different movement patterns:
#   - Directed: actively migrating cells (high confinement ratio, low turning)
#   - Confined: cells in dense tissue (low confinement ratio)
#   - Random: cells exploring their environment
#
#   We classify based on confinement ratio and linearity.
#

cat("\n=== PART 9: MOVEMENT TYPE CLASSIFICATION ===\n")

# Classify tracks based on movement metrics
tracks_classified <- tracks_clean %>%
  mutate(
    movement_type = case_when(
      CONFINEMENT_RATIO > 0.6 & LINEARITY_OF_FORWARD_PROGRESSION > 0.5 ~ "Directed",
      CONFINEMENT_RATIO < 0.3 ~ "Confined",
      TRUE ~ "Random"
    ),
    movement_type = factor(movement_type, levels = c("Directed", "Random", "Confined"))
  )

# Summary
movement_summary <- tracks_classified %>%
  count(movement_type) %>%
  mutate(pct = 100 * n / sum(n))

cat("Movement Type Classification:\n")
for (i in 1:nrow(movement_summary)) {
  cat(sprintf("  %s: %d tracks (%.1f%%)\n",
              movement_summary$movement_type[i],
              movement_summary$n[i],
              movement_summary$pct[i]))
}

# Plot 12: Scatter of confinement vs linearity
#
# INTERPRETATION:
#   - Upper right = directed migration
#   - Lower left = confined/exploratory
#   - Middle = random/intermediate
#
p12 <- ggplot(tracks_classified,
              aes(x = CONFINEMENT_RATIO, y = LINEARITY_OF_FORWARD_PROGRESSION,
                  color = movement_type)) +
  geom_point(data = tracks_classified %>% filter(movement_type == "Random"),
             alpha = 0.1, size = 0.5) +
  geom_point(data = tracks_classified %>% filter(movement_type != "Random"),
             alpha = 0.5, size = 1) +
  geom_vline(xintercept = c(0.3, 0.6), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = movement_colors, breaks = c("Directed", "Random", "Confined")) +
  labs(
    title = "Movement Type Classification",
    subtitle = "Based on confinement ratio and linearity",
    x = "Confinement ratio (displacement/distance)",
    y = "Linearity of forward progression",
    color = NULL
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
  theme_track()

# Plot 13: Movement type proportions (bar plot)
p13 <- ggplot(movement_summary, aes(x = movement_type, y = pct, fill = movement_type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
            vjust = -0.2, size = 3, fontface = "bold") +
  scale_fill_manual(values = movement_colors, breaks = c("Directed", "Random", "Confined")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Movement Type Distribution",
    x = NULL,
    y = "% of tracks",
    fill = NULL
  ) +
  theme_track() +
  theme(legend.position = "none", panel.grid.major.x = element_blank())

# Combine classification plots
classification_combined <- (p12 + p13) +
  plot_annotation(
    title = "Nuclear Movement Classification",
    subtitle = "Directed = persistent migration; Confined = restricted movement; Random = exploratory",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
    )
  )

ggsave("analysis_output/06_movement_classification.pdf", classification_combined, width = 12, height = 5)
cat("  Saved: 06_movement_classification.pdf\n")

# =============================================================================
# PART 10: CLUSTERING ANALYSIS (memory-efficient)
# =============================================================================
#
# BIOLOGICAL INTERPRETATION:
#   Clustering identifies subpopulations of nuclei with similar movement:
#   - Different cell types may have distinct motility signatures
#   - Regional differences in tissue dynamics
#   - Developmental stage-specific behaviors
#
#   We use k-means on a sample to avoid memory issues.
#

cat("\n=== PART 10: CLUSTERING ANALYSIS ===\n")

# Prepare features for clustering
clustering_features <- tracks_classified %>%
  select(TRACK_ID, SPEED_MEAN_UM_MIN, CONFINEMENT_RATIO,
         LINEARITY_OF_FORWARD_PROGRESSION, MEAN_DIRECTIONAL_CHANGE_RATE,
         DURATION_MIN) %>%
  filter(complete.cases(.))

# Sample for memory efficiency
set.seed(42)
sample_size <- min(8000, nrow(clustering_features))
sample_idx <- sample(1:nrow(clustering_features), sample_size)
features_sample <- clustering_features[sample_idx, ]

# Scale features
features_scaled <- scale(features_sample %>% select(-TRACK_ID))

# Determine optimal k using within-cluster sum of squares (faster than silhouette)
wss <- sapply(2:8, function(k) {
  km <- kmeans(features_scaled, centers = k, nstart = 10, iter.max = 20)
  km$tot.withinss
})

# Elbow method - look for bend
optimal_k <- 4  # Default to 4 clusters

cat(sprintf("Using %d clusters based on elbow method\n", optimal_k))

# Final clustering
set.seed(42)
km_result <- kmeans(features_scaled, centers = optimal_k, nstart = 25)
features_sample$cluster <- factor(km_result$cluster)

# Cluster summary
cluster_summary <- features_sample %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    mean_speed = mean(SPEED_MEAN_UM_MIN),
    mean_confinement = mean(CONFINEMENT_RATIO),
    mean_linearity = mean(LINEARITY_OF_FORWARD_PROGRESSION),
    mean_duration = mean(DURATION_MIN)
  ) %>%
  mutate(pct = 100 * n / sum(n))

cat("\nCluster Summary:\n")
print(as.data.frame(cluster_summary))

# Plot 14: Cluster characteristics
cluster_long <- features_sample %>%
  select(cluster, SPEED_MEAN_UM_MIN, CONFINEMENT_RATIO, LINEARITY_OF_FORWARD_PROGRESSION) %>%
  pivot_longer(-cluster, names_to = "metric") %>%
  mutate(metric = recode(metric,
                         "SPEED_MEAN_UM_MIN" = "Speed (µm/min)",
                         "CONFINEMENT_RATIO" = "Confinement ratio",
                         "LINEARITY_OF_FORWARD_PROGRESSION" = "Linearity"))

p14 <- ggplot(cluster_long, aes(x = cluster, y = value, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_viridis_d() +
  labs(
    title = "Cluster Characteristics",
    subtitle = "Movement properties by cluster",
    x = "Cluster",
    y = "Value"
  ) +
  theme_track() +
  theme(legend.position = "none")

ggsave("analysis_output/07_clustering.pdf", p14, width = 10, height = 5)
cat("  Saved: 07_clustering.pdf\n")

# =============================================================================
# PART 11: CONFINEMENT AND DISPLACEMENT ANALYSIS
# =============================================================================
#
# BIOLOGICAL INTERPRETATION:
#   Confinement ratio = net displacement / total distance traveled
#   - CR = 1: perfectly straight path (ballistic)
#   - CR → 0: highly confined (returns to origin)
#
#   In embryos:
#   - Migrating cells: CR > 0.5
#   - Cells in dense tissue: CR < 0.3
#   - Random exploration: CR ~ 0.3-0.5
#

cat("\n=== PART 11: CONFINEMENT ANALYSIS ===\n")

# Plot 15: Confinement ratio distribution
#
# INTERPRETATION:
#   - Left peak: confined nuclei (surrounded by neighbors)
#   - Right peak: migrating nuclei (directed movement)
#   - Spread: heterogeneity in cell behaviors
#
p15 <- ggplot(tracks_clean, aes(x = CONFINEMENT_RATIO)) +
  geom_histogram(bins = 50, fill = "#2166AC", alpha = 0.7) +
  geom_vline(xintercept = c(0.3, 0.7), linetype = "dashed", color = "#B2182B", linewidth = 0.4) +
  labs(
    title = "Confinement Ratio Distribution",
    subtitle = "<0.3 = confined; >0.7 = directed; between = random",
    x = "Confinement ratio",
    y = "Count"
  ) +
  theme_track()

# Plot 16: Confinement vs speed
#
# INTERPRETATION:
#   - Fast + high confinement = directed migration
#   - Slow + low confinement = stationary/confined
#   - Correlation indicates coupling between speed and directionality
#
p16 <- ggplot(tracks_clean, aes(x = CONFINEMENT_RATIO, y = SPEED_MEAN_UM_MIN)) +
  geom_point(alpha = 0.2, size = 0.5, color = "#2166AC") +
  geom_smooth(method = "loess", color = "#B2182B", linewidth = 1, se = FALSE) +
  labs(
    title = "Confinement vs Speed",
    subtitle = "Are faster nuclei more directed?",
    x = "Confinement ratio",
    y = "Mean speed (µm/min)"
  ) +
  theme_track()

# Combine confinement plots
confinement_combined <- (p15 + p16) +
  plot_annotation(
    title = "Confinement Analysis",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

ggsave("analysis_output/08_confinement.pdf", confinement_combined, width = 12, height = 5)
cat("  Saved: 08_confinement.pdf\n")

# =============================================================================
# PART 12: Z-DEPTH ANALYSIS
# =============================================================================
#
# BIOLOGICAL INTERPRETATION:
#   Movement may vary with depth in the embryo:
#   - Surface cells may migrate differently than deep cells
#   - Tissue boundaries affect movement
#   - Different cell types stratify by depth
#

cat("\n=== PART 12: Z-DEPTH ANALYSIS ===\n")

# Bin by Z-depth
z_analysis <- spots_velocity %>%
  filter(!is.na(inst_speed_um_min)) %>%
  mutate(z_layer = cut(POSITION_Z, breaks = 5,
                       labels = c("Deep", "Mid-deep", "Middle", "Mid-surface", "Surface")))

z_summary <- z_analysis %>%
  group_by(z_layer) %>%
  summarise(
    mean_speed = mean(inst_speed_um_min, na.rm = TRUE),
    se_speed = sd(inst_speed_um_min, na.rm = TRUE) / sqrt(n()),
    n = n()
  )

# Plot 17: Speed by Z-layer
#
# INTERPRETATION:
#   - Differences reveal tissue organization
#   - Surface vs deep movement patterns
#   - May reflect different cell populations
#
p17 <- ggplot(z_analysis, aes(x = z_layer, y = inst_speed_um_min, fill = z_layer)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.1, outlier.size = 0.3) +
  scale_fill_viridis_d() +
  coord_cartesian(ylim = c(0, 10)) +
  labs(
    title = "Nuclear Speed by Z-Depth (3D speed)",
    subtitle = "3D speed at different depths",
    x = "Z-layer (deep to surface)",
    y = "Instantaneous speed (µm/min)"
  ) +
  theme_track() +
  theme(legend.position = "none")

ggsave("analysis_output/09_z_depth_basic.pdf", p17, width = 8, height = 5)
cat("  Saved: 09_z_depth_basic.pdf\n")

# =============================================================================
# PART 12B: COMPREHENSIVE 4D SPATIAL-TEMPORAL ANALYSIS
# =============================================================================
#
# BIOLOGICAL INTERPRETATION:
#   Embryonic development is inherently 4D (X, Y, Z, Time). Nuclear movement
#   patterns change across:
#   - SPACE: Different regions (anterior/posterior, dorsal/ventral, deep/surface)
#   - TIME: Different developmental stages (early divisions vs later morphogenesis)
#   - COUPLED SPACE-TIME: Waves of movement, morphogenetic flows
#
#   This comprehensive analysis explores:
#   1. Z-dependent movement in XY space
#   2. Temporal evolution of spatial patterns
#   3. Full 4D parameter exploration
#   4. Interactive/3D visualizations
#

cat("\n=== PART 12B: COMPREHENSIVE 4D ANALYSIS ===\n")

# -----------------------------------------------------------------------------
# 12B.1: Spatial binning for 4D analysis
# -----------------------------------------------------------------------------

# Create spatial bins (voxels) for local statistics
SPATIAL_BIN_SIZE <- 15  # µm per bin
TIME_BIN_SIZE <- 10     # minutes per time bin

# Join velocity data with turning angle from spots_direction
spots_4d <- spots_velocity %>%
  filter(!is.na(inst_speed_um_min)) %>%
  # Add turning angle from spots_direction (calculated in Part 6)
  left_join(
    spots_direction %>%
      select(TRACK_ID, FRAME, turning_angle, direction_xy) %>%
      filter(!is.na(turning_angle)),
    by = c("TRACK_ID", "FRAME")
  ) %>%
  mutate(
    # Spatial bins
    x_bin = floor(POSITION_X / SPATIAL_BIN_SIZE) * SPATIAL_BIN_SIZE,
    y_bin = floor(POSITION_Y / SPATIAL_BIN_SIZE) * SPATIAL_BIN_SIZE,
    z_bin = floor(POSITION_Z / SPATIAL_BIN_SIZE) * SPATIAL_BIN_SIZE,

    # Time bins (in minutes)
    time_min = FRAME * FRAME_INTERVAL_MIN,
    time_bin = floor(time_min / TIME_BIN_SIZE) * TIME_BIN_SIZE,

    # Z-layers (normalized to 0-1 scale)
    z_range = max(POSITION_Z, na.rm = TRUE) - min(POSITION_Z, na.rm = TRUE),
    z_normalized = (POSITION_Z - min(POSITION_Z, na.rm = TRUE)) / z_range,
    z_layer = cut(z_normalized, breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                  labels = c("Deep (0-25%)", "Mid-deep (25-50%)",
                             "Mid-surface (50-75%)",
                            "Surface (75-100%)"),
                  include.lowest = TRUE)
  )

# -----------------------------------------------------------------------------
# 12B.2: Z-layer specific XY velocity patterns
# -----------------------------------------------------------------------------
#
# INTERPRETATION:
#   Different Z-layers may show different flow patterns:
#   - Deep cells: may be more stationary (yolk syncytial layer)
#   - Middle cells: active gastrulation movements
#   - Surface cells: epiboly, directed migration
#

cat("  Computing Z-layer specific XY patterns...\n")

# Local velocity vectors for each Z-layer
velocity_by_z <- spots_4d %>%
  filter(!is.na(vx) & !is.na(vy) & !is.na(z_layer)) %>%
  group_by(z_layer, x_bin, y_bin) %>%
  summarise(
    mean_speed = mean(inst_speed_um_min, na.rm = TRUE),
    mean_vx = mean(vx, na.rm = TRUE),
    mean_vy = mean(vy, na.rm = TRUE),
    mean_vz = mean(vz, na.rm = TRUE),
    sd_speed = sd(inst_speed_um_min, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 10)  # Need enough points for reliable statistics

# Plot: XY velocity heatmaps for each Z-layer
p_z_xy_heat <- ggplot(velocity_by_z, aes(x = x_bin, y = y_bin, fill = mean_speed)) +
  geom_tile() +
  facet_wrap(~z_layer, ncol = 5) +
  scale_fill_viridis_c(name = "3D Speed\n(µm/min)", option = "plasma", limits = c(0, NA)) +
  coord_fixed() +
  labs(
    title = "XY Speed Distribution by Z-Layer (3D speed)",
    subtitle = "3D speed mapped to XY positions at each depth",
    x = "X position (µm)",
    y = "Y position (µm)"
  ) +
  theme_track() +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "right"
  )

# Plot: XY velocity vectors for each Z-layer
p_z_xy_vector <- ggplot(velocity_by_z %>% filter(n >= 20),
                        aes(x = x_bin, y = y_bin)) +
  geom_segment(aes(xend = x_bin + mean_vx * 5, yend = y_bin + mean_vy * 5,
                   color = mean_speed),
               arrow = arrow(length = unit(0.1, "cm")), linewidth = 0.3) +
  facet_wrap(~z_layer, ncol = 5) +
  scale_color_viridis_c(name = "3D Speed\n(µm/min)", option = "plasma") +
  coord_fixed() +
  labs(
    title = "XY Velocity Vectors by Z-Layer",
    subtitle = "XY flow directions (Z-component not shown); arrows scaled 5x",
    x = "X position (µm)",
    y = "Y position (µm)"
  ) +
  theme_track() +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "right"
  )

# Combine Z-layer XY plots
z_xy_combined <- p_z_xy_heat / p_z_xy_vector +
  plot_annotation(
    title = "Depth-Dependent Movement Patterns in XY",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

ggsave("analysis_output/10_z_layer_xy_patterns.pdf", z_xy_combined, width = 14, height = 12)
cat("  Saved: 10_z_layer_xy_patterns.pdf\n")

# -----------------------------------------------------------------------------
# 12B.3: Z-movement analysis (vertical migration)
# -----------------------------------------------------------------------------
#
# INTERPRETATION:
#   Vertical (Z) movement is often overlooked but critical:
#   - Interkinetic nuclear migration: nuclei move up/down in neuroepithelium
#   - Involution during gastrulation: cells move from surface to deep
#   - Cell division: nuclei move apically before dividing
#

# Analyze vertical velocity component
z_movement_stats <- spots_4d %>%
  filter(!is.na(vz)) %>%
  group_by(z_layer) %>%
  summarise(
    mean_vz = mean(vz, na.rm = TRUE),
    se_vz = sd(vz, na.rm = TRUE) / sqrt(n()),
    prop_moving_up = mean(vz > 0, na.rm = TRUE),
    prop_moving_down = mean(vz < 0, na.rm = TRUE),
    n = n()
  )

# Plot: Z-velocity by layer
p_vz_layer <- ggplot(z_movement_stats, aes(x = z_layer, y = mean_vz, fill = mean_vz > 0)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_vz - se_vz, ymax = mean_vz + se_vz), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                    labels = c("TRUE" = "Moving up", "FALSE" = "Moving down")) +
  labs(
    title = "Vertical (Z-only) Velocity by Depth",
    subtitle = "+Z = towards surface; -Z = towards deep (Z-component only)",
    x = "Z-layer",
    y = "Mean Z-velocity (µm/min)",
    fill = NULL
  ) +
  theme_track() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Z-velocity distribution at each layer
p_vz_dist <- ggplot(spots_4d %>% filter(!is.na(vz) & abs(vz) < 10),
                    aes(x = vz, fill = z_layer)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~z_layer, ncol = 1, scales = "free_y") +
  scale_fill_viridis_d() +
  labs(
    title = "Z-Velocity Distribution by Layer",
    subtitle = "Are nuclei at certain depths preferentially moving up or down?",
    x = "Z-velocity (µm/min)",
    y = "Density"
  ) +
  theme_track() +
  theme(legend.position = "none",
        strip.text = element_text(size = 8))

z_movement_combined <- (p_vz_layer | p_vz_dist) +
  plot_annotation(
    title = "Vertical Nuclear Migration Analysis",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave("analysis_output/11_z_velocity_analysis.pdf", z_movement_combined, width = 12, height = 8)
cat("  Saved: 11_z_velocity_analysis.pdf\n")

# -----------------------------------------------------------------------------
# 12B.4: Temporal evolution of spatial patterns
# -----------------------------------------------------------------------------
#
# INTERPRETATION:
#   Movement patterns evolve during development:
#   - Early: more uniform, exploratory movement
#   - Later: coordinated morphogenetic flows
#   - Waves of activity may propagate across the embryo
#

cat("  Computing temporal evolution...\n")

# Speed over time for each Z-layer
speed_time_z <- spots_4d %>%
  filter(!is.na(inst_speed_um_min)) %>%
  group_by(time_bin, z_layer) %>%
  summarise(
    mean_speed = mean(inst_speed_um_min, na.rm = TRUE),
    se_speed = sd(inst_speed_um_min, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 20)

# Plot: Speed vs time for each Z-layer
p_speed_time_z <- ggplot(speed_time_z, aes(x = time_bin, y = mean_speed,
                                            color = z_layer, fill = z_layer)) +
  geom_ribbon(aes(ymin = mean_speed - se_speed, ymax = mean_speed + se_speed),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(
    title = "Nuclear Speed Over Time by Z-Layer",
    subtitle = "How movement activity evolves at different depths",
    x = "Time (min)",
    y = "Mean speed (µm/min)",
    color = "Z-layer",
    fill = "Z-layer"
  ) +
  theme_track() +
  theme(legend.position = "right")

# Speed evolution in XY space over time (sampled timepoints)
time_points <- unique(spots_4d$time_bin)
sample_times <- time_points[seq(1, length(time_points), length.out = min(6, length(time_points)))]

speed_xy_time <- spots_4d %>%
  filter(time_bin %in% sample_times) %>%
  group_by(time_bin, x_bin, y_bin) %>%
  summarise(
    mean_speed = mean(inst_speed_um_min, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 5)

p_speed_xy_time <- ggplot(speed_xy_time, aes(x = x_bin, y = y_bin, fill = mean_speed)) +
  geom_tile() +
  facet_wrap(~paste0("t = ", time_bin, " min"), ncol = 6) +
  scale_fill_viridis_c(name = "Speed\n(µm/min)", option = "plasma", limits = c(0, NA)) +
  coord_fixed() +
  labs(
    title = "Spatial Speed Patterns Over Time",
    subtitle = "How movement hotspots shift during development",
    x = "X position (µm)",
    y = "Y position (µm)"
  ) +
  theme_track() +
  theme(strip.text = element_text(face = "bold"))

# Combine temporal plots
temporal_combined <- p_speed_time_z / p_speed_xy_time +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(
    title = "Temporal Evolution of Movement Patterns",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave("analysis_output/12_temporal_evolution.pdf", temporal_combined, width = 12, height = 14)
cat("  Saved: 12_temporal_evolution.pdf\n")

# -----------------------------------------------------------------------------
# 12B.5: XZ and YZ projections (sagittal and coronal views)
# -----------------------------------------------------------------------------
#
# INTERPRETATION:
#   XY projection alone misses depth information.
#   XZ (sagittal) and YZ (coronal) views reveal:
#   - Vertical tissue organization
#   - Morphogenetic flows in the Z-axis
#   - Layer-specific movement patterns
#

# XZ projection
speed_xz <- spots_4d %>%
  group_by(x_bin, z_bin) %>%
  summarise(
    mean_speed = mean(inst_speed_um_min, na.rm = TRUE),
    mean_vx = mean(vx, na.rm = TRUE),
    mean_vz = mean(vz, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 10)

# YZ projection
speed_yz <- spots_4d %>%
  group_by(y_bin, z_bin) %>%
  summarise(
    mean_speed = mean(inst_speed_um_min, na.rm = TRUE),
    mean_vy = mean(vy, na.rm = TRUE),
    mean_vz = mean(vz, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 10)

# XZ heatmap and vectors
p_xz_heat <- ggplot(speed_xz, aes(x = x_bin, y = z_bin, fill = mean_speed)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Speed\n(µm/min)", option = "plasma") +
  labs(
    title = "XZ Projection (Sagittal View)",
    subtitle = "Speed distribution in anterior-posterior × depth",
    x = "X position (µm)",
    y = "Z position (µm)"
  ) +
  coord_fixed() +
  theme_track()

p_xz_vector <- ggplot(speed_xz %>% filter(n >= 20),
                      aes(x = x_bin, y = z_bin)) +
  geom_segment(aes(xend = x_bin + mean_vx * 5, yend = z_bin + mean_vz * 5,
                   color = mean_speed),
               arrow = arrow(length = unit(0.1, "cm")), linewidth = 0.4) +
  scale_color_viridis_c(name = "Speed\n(µm/min)", option = "plasma") +
  labs(
    title = "XZ Velocity Vectors",
    subtitle = "Flow patterns in sagittal plane",
    x = "X position (µm)",
    y = "Z position (µm)"
  ) +
  coord_fixed() +
  theme_track()

# YZ heatmap and vectors
p_yz_heat <- ggplot(speed_yz, aes(x = y_bin, y = z_bin, fill = mean_speed)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Speed\n(µm/min)", option = "plasma") +
  labs(
    title = "YZ Projection (Coronal View)",
    subtitle = "Speed distribution in left-right × depth",
    x = "Y position (µm)",
    y = "Z position (µm)"
  ) +
  coord_fixed() +
  theme_track()

p_yz_vector <- ggplot(speed_yz %>% filter(n >= 20),
                      aes(x = y_bin, y = z_bin)) +
  geom_segment(aes(xend = y_bin + mean_vy * 5, yend = z_bin + mean_vz * 5,
                   color = mean_speed),
               arrow = arrow(length = unit(0.1, "cm")), linewidth = 0.4) +
  scale_color_viridis_c(name = "Speed\n(µm/min)", option = "plasma") +
  labs(
    title = "YZ Velocity Vectors",
    subtitle = "Flow patterns in coronal plane",
    x = "Y position (µm)",
    y = "Z position (µm)"
  ) +
  coord_fixed() +
  theme_track()

# Combine projections
projections_combined <- (p_xz_heat + p_xz_vector) / (p_yz_heat + p_yz_vector) +
  plot_annotation(
    title = "Orthogonal Projections: XZ and YZ Views",
    subtitle = "Side and front views of 3D movement patterns",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey50")
    )
  )

ggsave("analysis_output/13_xz_yz_projections.pdf", projections_combined, width = 14, height = 12)
cat("  Saved: 13_xz_yz_projections.pdf\n")

# -----------------------------------------------------------------------------
# 12B.6: 4D Parameter Exploration - Correlations
# -----------------------------------------------------------------------------
#
# INTERPRETATION:
#   Understanding which parameters co-vary helps identify:
#   - Movement signatures for different cell populations
#   - Coupling between speed, directionality, and position
#   - Developmental gradients
#

cat("  Computing 4D parameter correlations...\n")

# Per-track summary with spatial info
track_spatial_summary <- spots_4d %>%
  group_by(TRACK_ID) %>%
  summarise(
    # Position (mean)
    mean_x = mean(POSITION_X, na.rm = TRUE),
    mean_y = mean(POSITION_Y, na.rm = TRUE),
    mean_z = mean(POSITION_Z, na.rm = TRUE),
    z_normalized = mean(z_normalized, na.rm = TRUE),

    # Time
    start_time = min(time_min, na.rm = TRUE),

    # Movement
    mean_speed = mean(inst_speed_um_min, na.rm = TRUE),
    mean_vz = mean(vz, na.rm = TRUE),
    speed_variability = sd(inst_speed_um_min, na.rm = TRUE) / mean(inst_speed_um_min, na.rm = TRUE),

    # Direction
    mean_turning = mean(turning_angle, na.rm = TRUE),

    n_spots = n(),
    .groups = "drop"
  ) %>%
  filter(n_spots >= 5) %>%
  left_join(
    tracks_classified %>% select(TRACK_ID, CONFINEMENT_RATIO, LINEARITY_OF_FORWARD_PROGRESSION,
                            TRACK_DISPLACEMENT, movement_type),
    by = "TRACK_ID"
  )

# Correlation matrix
cor_vars <- track_spatial_summary %>%
  select(mean_x, mean_y, mean_z, start_time, mean_speed, mean_vz,
         speed_variability, CONFINEMENT_RATIO, LINEARITY_OF_FORWARD_PROGRESSION)

cor_matrix <- cor(cor_vars, use = "pairwise.complete.obs")

# Plot correlation heatmap
cor_long <- as.data.frame(cor_matrix) %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "correlation") %>%
  mutate(
    var1 = factor(var1, levels = colnames(cor_matrix)),
    var2 = factor(var2, levels = rev(colnames(cor_matrix)))
  )

p_cor <- ggplot(cor_long, aes(x = var1, y = var2, fill = correlation)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  labs(
    title = "Parameter Correlation Matrix",
    subtitle = "Relationships between spatial, temporal, and movement features",
    x = NULL, y = NULL
  ) +
  theme_track() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

ggsave("analysis_output/14_parameter_correlations.pdf", p_cor, width = 10, height = 8)
cat("  Saved: 14_parameter_correlations.pdf\n")

# -----------------------------------------------------------------------------
# 12B.7: Movement type spatial distribution
# -----------------------------------------------------------------------------
#
# INTERPRETATION:
#   Where in the embryo do we find directed vs confined vs random movement?
#   - Directed: actively migrating regions
#   - Confined: densely packed regions
#   - Random: exploratory behavior
#

# Spatial distribution of movement types
movement_spatial <- track_spatial_summary %>%
  filter(!is.na(movement_type))

# XY distribution by movement type
p_movetype_xy <- ggplot(movement_spatial, aes(x = mean_x, y = mean_y, color = movement_type)) +
  geom_point(alpha = 0.3, size = 0.5) +
  facet_wrap(~movement_type) +
  scale_color_manual(values = movement_colors) +
  coord_fixed() +
  labs(
    title = "Spatial Distribution of Movement Types (XY)",
    subtitle = "Where in the embryo do we find each movement pattern?",
    x = "Mean X position (µm)",
    y = "Mean Y position (µm)"
  ) +
  theme_track() +
  theme(legend.position = "none")

# XZ distribution by movement type (showing Z-location for each type)
p_movetype_xz <- ggplot(movement_spatial, aes(x = mean_x, y = mean_z, color = movement_type)) +
  geom_point(alpha = 0.3, size = 0.5) +
  facet_wrap(~movement_type) +
  scale_color_manual(values = movement_colors) +
  labs(
    title = "Spatial Distribution of Movement Types (XZ - Side View)",
    subtitle = "Depth distribution of each movement pattern",
    x = "Mean X position (µm)",
    y = "Mean Z position (µm)"
  ) +
  theme_track() +
  theme(legend.position = "none")

# Z distribution by movement type (violin + boxplot)
p_movetype_z <- ggplot(movement_spatial, aes(x = movement_type, y = z_normalized,
                                              fill = movement_type)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.size = 0.3) +
  scale_fill_manual(values = movement_colors) +
  labs(
    title = "Z-Distribution by Movement Type",
    subtitle = "Are different movement types stratified by depth?",
    x = NULL,
    y = "Normalized Z position (0=deep, 1=surface)"
  ) +
  theme_track() +
  theme(legend.position = "none")

# Z density distribution (overlaid for comparison)
p_movetype_z_density <- ggplot(movement_spatial, aes(x = z_normalized, fill = movement_type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = movement_colors) +
  labs(
    title = "Z-Depth Density by Movement Type",
    subtitle = "Overlaid distributions show depth preferences",
    x = "Normalized Z position (0=deep, 1=surface)",
    y = "Density"
  ) +
  theme_track()

# Time distribution by movement type
p_movetype_time <- ggplot(movement_spatial, aes(x = start_time, fill = movement_type)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~movement_type, ncol = 1) +
  scale_fill_manual(values = movement_colors) +
  labs(
    title = "Temporal Distribution of Movement Types",
    subtitle = "When do different movement patterns occur?",
    x = "Track start time (min)",
    y = "Density"
  ) +
  theme_track() +
  theme(legend.position = "none")

# Summary statistics for Z by movement type
z_by_movetype_stats <- movement_spatial %>%
  group_by(movement_type) %>%
  summarise(
    n = n(),
    mean_z = mean(z_normalized, na.rm = TRUE),
    median_z = median(z_normalized, na.rm = TRUE),
    sd_z = sd(z_normalized, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n  Z-distribution by movement type:\n")
for (i in 1:nrow(z_by_movetype_stats)) {
  cat(sprintf("    %s: mean Z=%.2f, median=%.2f (n=%d)\n",
              z_by_movetype_stats$movement_type[i],
              z_by_movetype_stats$mean_z[i],
              z_by_movetype_stats$median_z[i],
              z_by_movetype_stats$n[i]))
}

movetype_spatial_combined <- (p_movetype_xy | p_movetype_xz) / (p_movetype_z | p_movetype_z_density) / p_movetype_time +
  plot_layout(heights = c(1.2, 1, 1.2)) +
  plot_annotation(
    title = "4D Distribution of Movement Types",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave("analysis_output/15_movement_type_4d.pdf", movetype_spatial_combined, width = 14, height = 16)
cat("  Saved: 15_movement_type_4d.pdf\n")

# -----------------------------------------------------------------------------
# 12B.8: 3D Track Visualization Data Export
# -----------------------------------------------------------------------------
#
# For interactive 3D visualization, we export data in formats compatible with:
#   - plotly (R/Python)
#   - napari (Python)
#   - ParaView (standalone)
#   - Fiji/ImageJ 3D viewer
#

cat("  Exporting 3D visualization data...\n")

# Sample tracks for 3D visualization (too many = cluttered)
set.seed(42)
sample_track_ids <- sample(unique(spots_clean$TRACK_ID), min(500, n_distinct(spots_clean$TRACK_ID)))

tracks_3d_export <- spots_clean %>%
  filter(TRACK_ID %in% sample_track_ids) %>%
  arrange(TRACK_ID, FRAME) %>%
  left_join(
    tracks_classified %>% dplyr::select(TRACK_ID, SPEED_MEAN_UM_MIN, CONFINEMENT_RATIO, movement_type),
    by = "TRACK_ID"
  ) %>%
  select(TRACK_ID, FRAME, POSITION_X, POSITION_Y, POSITION_Z,
         SPEED_MEAN_UM_MIN, CONFINEMENT_RATIO, movement_type)

write_csv(tracks_3d_export, "analysis_output/tracks_3d_visualization.csv")
cat("    Exported: tracks_3d_visualization.csv (", nrow(tracks_3d_export), " points, ",
    n_distinct(tracks_3d_export$TRACK_ID), " tracks)\n")

# Also create a summary for plotly-style 3D scatter
track_centroids <- track_spatial_summary %>%
  filter(!is.na(movement_type)) %>%
  select(TRACK_ID, mean_x, mean_y, mean_z, mean_speed, CONFINEMENT_RATIO, movement_type)

write_csv(track_centroids, "analysis_output/track_centroids_3d.csv")

# -----------------------------------------------------------------------------
# 12B.9: Static 3D-like visualization (pseudo-3D)
# -----------------------------------------------------------------------------
#
# Create a multi-panel view that gives 3D impression:
#   - XY colored by Z
#   - XZ side view
#   - YZ front view
#   - Time-color animation frames
#

# Use ALL tracks for complete visualization (not sampled)
cat("  Creating pseudo-3D visualization with all tracks...\n")

viz_all <- spots_clean %>%
  arrange(TRACK_ID, FRAME) %>%
  mutate(
    z_normalized = (POSITION_Z - min(POSITION_Z)) / (max(POSITION_Z) - min(POSITION_Z)),
    time_normalized = FRAME / max(FRAME)
  )

# XY view colored by Z
p_3d_xy <- ggplot(viz_all, aes(x = POSITION_X, y = POSITION_Y,
                                group = TRACK_ID, color = z_normalized)) +
  geom_path(alpha = 0.3, linewidth = 0.15) +
  scale_color_viridis_c(name = "Z depth\n(normalized)", option = "plasma") +
  coord_fixed() +
  labs(title = "XY View (top-down)", subtitle = "Color = depth", x = "X (µm)", y = "Y (µm)") +
  theme_track()

# XZ view (sagittal)
p_3d_xz <- ggplot(viz_all, aes(x = POSITION_X, y = POSITION_Z,
                                group = TRACK_ID, color = time_normalized)) +
  geom_path(alpha = 0.3, linewidth = 0.15) +
  scale_color_viridis_c(name = "Time", option = "viridis", guide = "none") +
  labs(title = "XZ View (side)", subtitle = "Color = time", x = "X (µm)", y = "Z (µm)") +
  theme_track()

# YZ view (coronal)
p_3d_yz <- ggplot(viz_all, aes(x = POSITION_Y, y = POSITION_Z,
                                group = TRACK_ID, color = time_normalized)) +
  geom_path(alpha = 0.3, linewidth = 0.15) +
  scale_color_viridis_c(name = "Time", option = "viridis", guide = "none") +
  labs(title = "YZ View (front)", subtitle = "Color = time", x = "Y (µm)", y = "Z (µm)") +
  theme_track()

# XY view colored by time
p_3d_xy_time <- ggplot(viz_all, aes(x = POSITION_X, y = POSITION_Y,
                                     group = TRACK_ID, color = time_normalized)) +
  geom_path(alpha = 0.3, linewidth = 0.15) +
  scale_color_viridis_c(name = "Time\n(normalized)", option = "viridis") +
  coord_fixed() +
  labs(title = "XY View (top-down)", subtitle = "Color = time progression", x = "X (µm)", y = "Y (µm)") +
  theme_track()

# Combine into comprehensive pseudo-3D view
pseudo_3d <- (p_3d_xy | p_3d_xy_time) / (p_3d_xz | p_3d_yz) +
  plot_annotation(
    title = "Complete 3D Track Visualization - All Tracks",
    subtitle = sprintf("%d tracks, %d spots | Multiple viewing angles",
                      n_distinct(viz_all$TRACK_ID), nrow(viz_all)),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey50")
    )
  )

ggsave("analysis_output/16_pseudo_3d_tracks.pdf", pseudo_3d, width = 16, height = 14)
cat("  Saved: 16_pseudo_3d_tracks.pdf\n")

# =============================================================================
# PART 14: INTERACTIVE 4D VISUALIZATION WITH PLOTLY
# =============================================================================
#
# BIOLOGICAL INTERPRETATION:
#   Interactive 3D visualization allows:
#   - Rotating the embryo to see tracks from any angle
#   - Identifying spatial clusters of similar movement
#   - Following individual tracks through time
#   - Zooming into regions of interest
#
#   We create multiple HTML files for different views.
#

cat("\n=== PART 14: INTERACTIVE 4D VISUALIZATION ===\n")

# -----------------------------------------------------------------------------
# 14.1: Interactive 3D track plot (static snapshot with rotation)
# -----------------------------------------------------------------------------

cat("  Creating interactive 3D visualizations...\n")

# Sample tracks for interactive view (too many = slow rendering)
set.seed(42)
n_sample_interactive <- min(1000, n_distinct(spots_clean$TRACK_ID))
sample_ids_interactive <- sample(unique(spots_clean$TRACK_ID), n_sample_interactive)

viz_interactive <- spots_clean %>%
  filter(TRACK_ID %in% sample_ids_interactive) %>%
  arrange(TRACK_ID, FRAME) %>%
  left_join(
    tracks_classified %>% select(TRACK_ID, movement_type, SPEED_MEAN_UM_MIN),
    by = "TRACK_ID"
  ) %>%
  mutate(
    time_min = FRAME * FRAME_INTERVAL_MIN,
    movement_type = as.character(movement_type)
  )

# 3D scatter plot of track positions colored by movement type
fig_3d_movetype <- plot_ly(
  viz_interactive,
  x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
  color = ~movement_type,
  colors = c("Directed" = "#2166AC", "Random" = "#762A83", "Confined" = "#B2182B"),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 1.5, opacity = 0.4),
  hoverinfo = "text",
  text = ~paste("Track:", TRACK_ID, "<br>Type:", movement_type,
                "<br>Speed:", round(SPEED_MEAN_UM_MIN, 2), "µm/min")
) %>%
  plotly::layout(
    title = list(
      text = "3D Nuclear Positions by Movement Type",
      font = list(size = 16)
    ),
    scene = list(
      xaxis = list(title = "X (µm)"),
      yaxis = list(title = "Y (µm)"),
      zaxis = list(title = "Z (µm)"),
      aspectmode = "data"
    ),
    legend = list(title = list(text = "Movement Type"))
  )

htmlwidgets::saveWidget(fig_3d_movetype, "analysis_output/17_interactive_3d_movetype.html",
                        selfcontained = TRUE)
cat("    Saved: 17_interactive_3d_movetype.html\n")

# 3D scatter colored by speed
fig_3d_speed <- plot_ly(
  viz_interactive,
  x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
  color = ~SPEED_MEAN_UM_MIN,
  colors = viridis::plasma(100),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 1.5, opacity = 0.5),
  hoverinfo = "text",
  text = ~paste("Track:", TRACK_ID, "<br>Speed:", round(SPEED_MEAN_UM_MIN, 2), "µm/min")
) %>%
  plotly::layout(
    title = list(
      text = "3D Nuclear Positions by Speed",
      font = list(size = 16)
    ),
    scene = list(
      xaxis = list(title = "X (µm)"),
      yaxis = list(title = "Y (µm)"),
      zaxis = list(title = "Z (µm)"),
      aspectmode = "data"
    )
  )

htmlwidgets::saveWidget(fig_3d_speed, "analysis_output/18_interactive_3d_speed.html",
                        selfcontained = TRUE)
cat("    Saved: 18_interactive_3d_speed.html\n")

# -----------------------------------------------------------------------------
# 14.2: 3D track lines (showing actual trajectories)
# -----------------------------------------------------------------------------

# Create line traces for each track
cat("  Creating 3D track line visualization...\n")

# Use fewer tracks for line visualization (more complex rendering)
n_line_tracks <- min(1000, n_distinct(spots_clean$TRACK_ID))
line_track_ids <- sample(unique(spots_clean$TRACK_ID), n_line_tracks)

viz_lines <- spots_clean %>%
  filter(TRACK_ID %in% line_track_ids) %>%
  arrange(TRACK_ID, FRAME) %>%
  left_join(
    tracks_classified %>% select(TRACK_ID, movement_type, SPEED_MEAN_UM_MIN),
    by = "TRACK_ID"
  )


# Create 3D line plot
fig_3d_lines <- plot_ly()

for (tid in unique(viz_lines$TRACK_ID)) {
  track_data <- viz_lines %>% filter(TRACK_ID == tid)
  mvtype <- unique(track_data$movement_type)[1]
  color <- switch(as.character(mvtype),
                  "Directed" = "#2166AC",
                  "Confined" = "#B2182B",
                  "#762A83")

  fig_3d_lines <- fig_3d_lines %>%
    add_trace(
      data = track_data,
      x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
      type = "scatter3d",
      mode = "lines",
      line = list(width = 2, color = color),
      opacity = 0.6,
      name = mvtype,
      legendgroup = mvtype,
      showlegend = (tid == line_track_ids[1] ||
                    !mvtype %in% sapply(line_track_ids[1:(which(line_track_ids == tid)-1)],
                                        function(x) unique(viz_lines$movement_type[viz_lines$TRACK_ID == x])[1])),
      hoverinfo = "text",
      text = paste("Track:", tid, "| Type:", mvtype)
    )
}

fig_3d_lines <- fig_3d_lines %>%
  plotly::layout(
    title = list(
      text = sprintf("3D Track Trajectories (n=%d sample)", n_line_tracks),
      font = list(size = 16)
    ),
    scene = list(
      xaxis = list(title = "X (µm)"),
      yaxis = list(title = "Y (µm)"),
      zaxis = list(title = "Z (µm)"),
      aspectmode = "data"
    )
  )

htmlwidgets::saveWidget(fig_3d_lines, "analysis_output/19_interactive_3d_tracks.html",
                        selfcontained = TRUE)
cat("    Saved: 19_interactive_3d_tracks.html\n")

# -----------------------------------------------------------------------------
# 14.3: 4D Animation (3D + Time slider)
# -----------------------------------------------------------------------------

cat("  Creating 4D time-lapse animation...\n")

# Bin time for animation frames
time_bins <- seq(0, max(viz_interactive$time_min), by = 5)  # 5-minute frames

viz_animated <- viz_interactive %>%
  mutate(time_bin = cut(time_min, breaks = c(time_bins, Inf),
                        labels = paste0(time_bins, "-", c(time_bins[-1], max(time_min) + 1), " min"),
                        include.lowest = TRUE))

fig_4d <- plot_ly(
  viz_animated,
  x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
  color = ~movement_type,
  colors = c("Directed" = "#2166AC", "Random" = "#762A83", "Confined" = "#B2182B"),
  frame = ~time_bin,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 2, opacity = 0.5)
) %>%
  plotly::layout(
    title = list(
      text = "4D Nuclear Movement (use slider for time)",
      font = list(size = 16)
    ),
    scene = list(
      xaxis = list(title = "X (µm)"),
      yaxis = list(title = "Y (µm)"),
      zaxis = list(title = "Z (µm)"),
      aspectmode = "data"
    )
  ) %>%
  animation_opts(
    frame = 500,  # ms per frame
    transition = 300,
    redraw = FALSE
  ) %>%
  animation_slider(
    currentvalue = list(prefix = "Time: ")
  )

htmlwidgets::saveWidget(fig_4d, "analysis_output/20_interactive_4d_timelapse.html",
                        selfcontained = TRUE)
cat("    Saved: 20_interactive_4d_timelapse.html\n")

# =============================================================================
# PART 15: ADDITIONAL ANALYSES
# =============================================================================
#
# Additional analyses that provide deeper biological insights:
#

cat("\n=== PART 15: ADDITIONAL ANALYSES ===\n")

# -----------------------------------------------------------------------------
# 15.1: Track persistence (how long do tracks maintain direction?)
# -----------------------------------------------------------------------------

cat("  Computing track persistence...\n")

# Calculate directional autocorrelation
calculate_persistence <- function(df) {
  if (nrow(df) < 5) return(NA)

  # Get displacement vectors
  displacements <- df %>%
    arrange(FRAME) %>%
    mutate(
      dx = POSITION_X - lag(POSITION_X),
      dy = POSITION_Y - lag(POSITION_Y),
      dz = POSITION_Z - lag(POSITION_Z)
    ) %>%
    filter(!is.na(dx))

  if (nrow(displacements) < 3) return(NA)

  # Compute autocorrelation of direction
  # (cosine of angle between consecutive displacements)
  autocorr <- displacements %>%
    mutate(
      mag = sqrt(dx^2 + dy^2 + dz^2),
      dx_norm = dx / (mag + 1e-10),
      dy_norm = dy / (mag + 1e-10),
      dz_norm = dz / (mag + 1e-10),
      dot_next = dx_norm * lead(dx_norm) + dy_norm * lead(dy_norm) + dz_norm * lead(dz_norm)
    )

  mean(autocorr$dot_next, na.rm = TRUE)
}

# Calculate persistence for each track
track_persistence <- spots_clean %>%
  group_by(TRACK_ID) %>%
  group_modify(~tibble(persistence = calculate_persistence(.x))) %>%
  ungroup() %>%
  filter(!is.na(persistence)) %>%
  left_join(tracks_classified %>% select(TRACK_ID, movement_type), by = "TRACK_ID")

# Plot persistence by movement type
p_persistence <- ggplot(track_persistence %>% filter(!is.na(movement_type)),
                        aes(x = movement_type, y = persistence, fill = movement_type)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = movement_colors) +
  labs(
    title = "Directional Persistence by Movement Type",
    subtitle = "1 = perfectly persistent; 0 = random; -1 = reversing",
    x = NULL,
    y = "Persistence (autocorrelation)"
  ) +
  theme_track() +
  theme(legend.position = "none")

ggsave("analysis_output/21_persistence.pdf", p_persistence, width = 8, height = 6)
cat("    Saved: 21_persistence.pdf\n")

# -----------------------------------------------------------------------------
# 15.2: Speed-confinement phase space
# -----------------------------------------------------------------------------
#
# INTERPRETATION:
#   A 2D phase space of speed vs confinement reveals distinct populations:
#   - Fast + directed: active migration
#   - Slow + confined: stationary/crowded cells
#   - Fast + confined: oscillating cells (e.g., interkinetic migration)
#   - Slow + directed: slow persistent crawling
#

p_phase_space <- ggplot(tracks_classified,
                        aes(x = SPEED_MEAN_UM_MIN, y = CONFINEMENT_RATIO, color = movement_type)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_density_2d(linewidth = 0.3, alpha = 0.7) +
  scale_color_manual(values = movement_colors) +
  labs(
    title = "Speed-Confinement Phase Space",
    subtitle = "Contours show density; movement type classification overlaid",
    x = "Mean speed (µm/min, 3D)",
    y = "Confinement ratio",
    color = NULL
  ) +
  theme_track()

ggsave("analysis_output/22_phase_space.pdf", p_phase_space, width = 10, height = 8)
cat("    Saved: 22_phase_space.pdf\n")

# -----------------------------------------------------------------------------
# 15.3: Local density vs speed analysis
# -----------------------------------------------------------------------------
#
# INTERPRETATION:
#   Does nuclear density affect movement speed?
#   - High density → more crowded → slower movement?
#   - Or high activity regions → both high density and speed?
#

cat("  Computing local density effects...\n")

# Calculate local density for each spot
spots_with_density <- spots_clean %>%
  mutate(
    x_bin = floor(POSITION_X / 100) * 100,
    y_bin = floor(POSITION_Y / 100) * 100,
    z_bin = floor(POSITION_Z / 50) * 50,
    time_bin = floor(FRAME / 10) * 10
  )

local_density <- spots_with_density %>%
  group_by(x_bin, y_bin, z_bin, time_bin) %>%
  summarise(local_density = n(), .groups = "drop")

spots_density <- spots_with_density %>%
  left_join(local_density, by = c("x_bin", "y_bin", "z_bin", "time_bin")) %>%
  left_join(
    spots_velocity %>% select(TRACK_ID, FRAME, inst_speed_um_min),
    by = c("TRACK_ID", "FRAME")
  ) %>%
  filter(!is.na(inst_speed_um_min) & !is.na(local_density))

# Plot density vs speed
density_speed_summary <- spots_density %>%
  mutate(density_bin = cut(local_density, breaks = quantile(local_density, probs = seq(0, 1, 0.1)),
                           include.lowest = TRUE)) %>%
  group_by(density_bin) %>%
  summarise(
    mean_speed = mean(inst_speed_um_min, na.rm = TRUE),
    se_speed = sd(inst_speed_um_min, na.rm = TRUE) / sqrt(n()),
    mean_density = mean(local_density),
    n = n()
  ) %>%
  filter(!is.na(density_bin))

p_density_speed <- ggplot(density_speed_summary, aes(x = mean_density, y = mean_speed)) +
  geom_point(aes(size = n), color = "#2166AC", alpha = 0.7) +
  geom_smooth(method = "loess", color = "#B2182B", se = TRUE, alpha = 0.2) +
  labs(
    title = "Local Density vs Nuclear Speed",
    subtitle = "Does crowding affect movement? (loess fit with 95% CI)",
    x = "Local nuclear density (nuclei per voxel)",
    y = "Mean instantaneous speed (µm/min, 3D)",
    size = "n obs"
  ) +
  theme_track()

ggsave("analysis_output/23_density_vs_speed.pdf", p_density_speed, width = 10, height = 7)
cat("    Saved: 23_density_vs_speed.pdf\n")

# -----------------------------------------------------------------------------
# 15.4: Neighbor correlation analysis
# -----------------------------------------------------------------------------
#
# INTERPRETATION:
#   Do nearby nuclei move together (correlated movement)?
#   - High correlation = collective behavior / tissue flow
#   - Low correlation = independent cell movement
#

cat("  Computing neighbor velocity correlation...\n")

# Sample frames for neighbor analysis (computationally expensive)
sample_frames <- sample(unique(spots_velocity$FRAME), min(20, n_distinct(spots_velocity$FRAME)))

neighbor_corr <- map_dfr(sample_frames, function(f) {
  frame_data <- spots_velocity %>%
    filter(FRAME == f & !is.na(vx) & !is.na(vy) & !is.na(vz))

  if (nrow(frame_data) < 50) return(NULL)

  # Sample spots to avoid O(n²) computation
  sample_spots <- sample_n(frame_data, min(200, nrow(frame_data)))

  # For each spot, find neighbors and compute velocity correlation
  correlations <- map_dfr(1:nrow(sample_spots), function(i) {
    spot <- sample_spots[i, ]

    # Find neighbors within 50 µm
    neighbors <- frame_data %>%
      mutate(
        dist = sqrt((POSITION_X - spot$POSITION_X)^2 +
                    (POSITION_Y - spot$POSITION_Y)^2 +
                    (POSITION_Z - spot$POSITION_Z)^2)
      ) %>%
      filter(dist > 0 & dist < 50)

    if (nrow(neighbors) == 0) return(NULL)

    # Compute velocity correlation with neighbors
    v_self <- c(spot$vx, spot$vy, spot$vz)
    v_self_norm <- v_self / (sqrt(sum(v_self^2)) + 1e-10)

    neighbor_corrs <- neighbors %>%
      rowwise() %>%
      mutate(
        v_neigh = list(c(vx, vy, vz)),
        v_neigh_norm = list(unlist(v_neigh) / (sqrt(sum(unlist(v_neigh)^2)) + 1e-10)),
        corr = sum(v_self_norm * unlist(v_neigh_norm))
      ) %>%
      ungroup()

    tibble(
      distance_bin = cut(neighbor_corrs$dist, breaks = c(0, 10, 20, 30, 40, 50)),
      correlation = neighbor_corrs$corr
    )
  })

  correlations
})

# Plot neighbor velocity correlation by distance
if (nrow(neighbor_corr) > 0) {
  neighbor_summary <- neighbor_corr %>%
    filter(!is.na(distance_bin)) %>%
    group_by(distance_bin) %>%
    summarise(
      mean_corr = mean(correlation, na.rm = TRUE),
      se_corr = sd(correlation, na.rm = TRUE) / sqrt(n()),
      n = n()
    )

  p_neighbor <- ggplot(neighbor_summary, aes(x = distance_bin, y = mean_corr)) +
    geom_col(fill = "#2166AC", alpha = 0.7) +
    geom_errorbar(aes(ymin = mean_corr - se_corr, ymax = mean_corr + se_corr), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = "Neighbor Velocity Correlation",
      subtitle = "Do nearby nuclei move in similar directions? (1 = same, 0 = random, -1 = opposite)",
      x = "Distance to neighbor (µm)",
      y = "Mean velocity correlation"
    ) +
    theme_track()

  ggsave("analysis_output/24_neighbor_correlation.pdf", p_neighbor, width = 10, height = 6)
  cat("    Saved: 24_neighbor_correlation.pdf\n")
}

# -----------------------------------------------------------------------------
# 12B.10: 4D Summary Statistics
# -----------------------------------------------------------------------------

cat("\n4D Analysis Summary:\n")

# Spatial extent
spatial_extent <- spots_clean %>%
  summarise(
    x_range = max(POSITION_X) - min(POSITION_X),
    y_range = max(POSITION_Y) - min(POSITION_Y),
    z_range = max(POSITION_Z) - min(POSITION_Z)
  )

cat(sprintf("  Spatial extent: X=%.0f µm, Y=%.0f µm, Z=%.0f µm\n",
            spatial_extent$x_range, spatial_extent$y_range, spatial_extent$z_range))

# Speed by Z-layer
z_speed_summary <- spots_4d %>%
  filter(!is.na(z_layer)) %>%
  group_by(z_layer) %>%
  summarise(mean_speed = mean(inst_speed_um_min, na.rm = TRUE)) %>%
  arrange(z_layer)

cat("  Speed by Z-layer:\n")
for (i in 1:nrow(z_speed_summary)) {
  cat(sprintf("    %s: %.2f µm/min\n",
              z_speed_summary$z_layer[i], z_speed_summary$mean_speed[i]))
}

# Movement type spatial bias
movetype_z_summary <- movement_spatial %>%
  group_by(movement_type) %>%
  summarise(
    mean_z_norm = mean(z_normalized, na.rm = TRUE),
    n = n()
  )

cat("  Movement type depth bias (0=deep, 1=surface):\n")
for (i in 1:nrow(movetype_z_summary)) {
  cat(sprintf("    %s: z=%.2f (n=%d)\n",
              movetype_z_summary$movement_type[i],
              movetype_z_summary$mean_z_norm[i],
              movetype_z_summary$n[i]))
}

# =============================================================================
# PART 13: FINAL SUMMARY
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("ANALYSIS SUMMARY\n")
cat(strrep("=", 70), "\n\n")

cat("DATASET:\n")
cat(sprintf("  Filtered tracks: %d (from %d original)\n", nrow(tracks_clean), nrow(tracks)))
cat(sprintf("  Total imaging time: %.1f min (%.1f hours)\n",
            max(spots$FRAME) * FRAME_INTERVAL_MIN,
            max(spots$FRAME) * FRAME_INTERVAL_MIN / 60))
cat(sprintf("  Temporal resolution: %d sec/frame\n", FRAME_INTERVAL_SEC))

cat("\nNUCLEAR MOVEMENT:\n")
cat(sprintf("  Median speed: %.2f µm/min\n", median(tracks_clean$SPEED_MEAN_UM_MIN)))
cat(sprintf("  Median track duration: %.1f min\n", median(tracks_clean$DURATION_MIN)))
cat(sprintf("  Median displacement: %.1f µm\n", median(tracks_clean$TRACK_DISPLACEMENT)))

cat("\nMOVEMENT TYPES:\n")
for (i in 1:nrow(movement_summary)) {
  cat(sprintf("  %s: %.1f%%\n", movement_summary$movement_type[i], movement_summary$pct[i]))
}

cat("\nDIFFUSION ANALYSIS:\n")
cat(sprintf("  MSD exponent (α): %.2f → %s\n", alpha,
            ifelse(alpha < 0.8, "subdiffusive (confined)",
                   ifelse(alpha > 1.2, "superdiffusive (directed)", "normal diffusion"))))

cat("\nOUTPUT FILES:\n")
cat("  Quality Control & Basic Analysis:\n")
cat("    analysis_output/01_quality_control.pdf\n")
cat("    analysis_output/02_velocity_analysis.pdf\n")
cat("    analysis_output/03_directionality.pdf\n")
cat("    analysis_output/04_spatial_velocity.pdf\n")
cat("    analysis_output/05_msd_analysis.pdf\n")
cat("    analysis_output/06_movement_classification.pdf\n")
cat("    analysis_output/07_clustering.pdf\n")
cat("    analysis_output/08_confinement.pdf\n")
cat("  4D Spatial-Temporal Analysis:\n")
cat("    analysis_output/09_z_depth_basic.pdf\n")
cat("    analysis_output/10_z_layer_xy_patterns.pdf\n")
cat("    analysis_output/11_z_velocity_analysis.pdf\n")
cat("    analysis_output/12_temporal_evolution.pdf\n")
cat("    analysis_output/13_xz_yz_projections.pdf\n")
cat("    analysis_output/14_parameter_correlations.pdf\n")
cat("    analysis_output/15_movement_type_4d.pdf\n")
cat("    analysis_output/16_pseudo_3d_tracks.pdf\n")
cat("  Interactive 4D Visualizations (open in browser):\n")
cat("    analysis_output/17_interactive_3d_movetype.html\n")
cat("    analysis_output/18_interactive_3d_speed.html\n")
cat("    analysis_output/19_interactive_3d_tracks.html\n")
cat("    analysis_output/20_interactive_4d_timelapse.html\n")
cat("  Additional Analyses:\n")
cat("    analysis_output/21_persistence.pdf\n")
cat("    analysis_output/22_phase_space.pdf\n")
cat("    analysis_output/23_density_vs_speed.pdf\n")
cat("    analysis_output/24_neighbor_correlation.pdf\n")
cat("  Data Exports:\n")

# Save processed data
write_csv(tracks_classified, "analysis_output/tracks_analyzed.csv")
write_csv(spots_velocity, "analysis_output/spots_with_velocity.csv")
cat("    analysis_output/tracks_analyzed.csv\n")
cat("    analysis_output/spots_with_velocity.csv\n")
cat("    analysis_output/tracks_3d_visualization.csv\n")
cat("    analysis_output/track_centroids_3d.csv\n")

cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n")
