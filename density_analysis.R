# =============================================================================
# Density Analysis — Geodesic Surface Density & Nuclei Counts
# =============================================================================
#
# PURPOSE:
#   Publication-quality density comparison between zebrafish and medaka.
#   Replicates the geodesic density approach from embryo_viewer.py in R:
#     θ,φ → unit-sphere XYZ → per-frame cKDTree (RANN) with adaptive radius
#
# PLOTS PRODUCED (all as multi-panel PDFs):
#   density_01_nuclei_count.pdf         — nuclei count per frame (global + ROI)
#   density_02_cumulative.pdf           — cumulative nuclei over time
#   density_03_rate_of_change.pdf       — Δ nuclei per frame (global + ROI)
#   density_04_geodesic_density.pdf     — mean geodesic surface density over time
#   density_05_combined_dashboard.pdf   — all-in-one summary panel
#
# INPUTS (per species):
#   - oriented_spots.csv   (FRAME, THETA_DEG, PHI_DEG, IN_ROI, TRACK_IN_ROI)
#   - sphere_params.csv    (fitted sphere radius)
#
# =============================================================================

source("renv/activate.R")

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(viridis)
library(RANN)

# =============================================================================
# PARAMETERS
# =============================================================================

# NOTE: Frame intervals differ between species (ZF=120s, MK=30s).
# We use the POSITION_T column from TrackMate which is already calibrated
# in seconds for each species, so no hardcoded interval is needed.

OUTPUT_DIR <- "analysis_output"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

SPECIES <- list(
  zebrafish = list(
    name     = "Zebrafish",
    short    = "ZF",
    color    = "#2166AC",
    fill     = "#4393C3",
    light    = "#D1E5F0",
    data_dir = "analysis_output_zebrafish_05112025"
  ),
  medaka = list(
    name     = "Medaka",
    short    = "MK",
    color    = "#B2182B",
    fill     = "#D6604D",
    light    = "#FDDBC7",
    data_dir = "output_medaka_28052025"
  )
)

# Geodesic density: k for adaptive radius
K_NN         <- 16
SUBSAMPLE_N  <- 5000

# =============================================================================
# THEME
# =============================================================================

species_colors <- c("Zebrafish" = "#2166AC", "Medaka" = "#B2182B")
species_fills  <- c("Zebrafish" = "#4393C3", "Medaka" = "#D6604D")

theme_density <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold",
                                      size = base_size + 2, margin = margin(b = 4)),
      plot.subtitle    = element_text(hjust = 0.5, size = base_size - 1,
                                      color = "grey45", margin = margin(b = 8)),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.text      = element_text(size = base_size),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey90"),
      strip.text       = element_text(face = "bold", size = base_size),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1),
      plot.margin      = margin(10, 12, 6, 8)
    )
}

save_pdf <- function(plot, filename, width = 14, height = 8) {
  path <- file.path(OUTPUT_DIR, filename)
  ggsave(path, plot, width = width, height = height, device = "pdf")
  cat(sprintf("  -> Saved %s\n", path))
}

# =============================================================================
# DATA LOADING
# =============================================================================

cat("\n", strrep("═", 70), "\n")
cat("  DENSITY ANALYSIS — Zebrafish vs. Medaka\n")
cat(strrep("═", 70), "\n\n")

all_spots <- list()

for (sp_key in names(SPECIES)) {
  sp <- SPECIES[[sp_key]]
  cat(sprintf("Loading %s from %s/ ... ", sp$name, sp$data_dir))

  spots_path <- file.path(sp$data_dir, "oriented_spots.csv")
  sphere_path <- file.path(sp$data_dir, "sphere_params.csv")

  dt <- fread(spots_path, select = c("FRAME", "POSITION_T", "THETA_DEG", "PHI_DEG",
                                      "IN_ROI", "TRACK_IN_ROI"))
  # Coerce logical
  dt[, IN_ROI := as.logical(IN_ROI)]
  dt[, TRACK_IN_ROI := as.logical(TRACK_IN_ROI)]
  dt[, species := sp$name]
  # POSITION_T is in seconds (calibrated by TrackMate per species)
  dt[, time_min := POSITION_T / 60]

  # Derive per-species frame interval from the data
  frame_t <- unique(dt[, .(FRAME, POSITION_T)])[order(FRAME)]
  if (nrow(frame_t) > 1) {
    SPECIES[[sp_key]]$frame_interval_sec <- frame_t$POSITION_T[2] - frame_t$POSITION_T[1]
  } else {
    SPECIES[[sp_key]]$frame_interval_sec <- NA_real_
  }

  # Read sphere radius
  sph <- fread(sphere_path)
  SPECIES[[sp_key]]$radius <- as.numeric(sph$radius[1])

  cat(sprintf("%s spots, %d frames, dt=%.0fs, R = %.1f um\n",
              format(nrow(dt), big.mark = ","), uniqueN(dt$FRAME),
              SPECIES[[sp_key]]$frame_interval_sec,
              SPECIES[[sp_key]]$radius))

  all_spots[[sp_key]] <- dt
}

spots <- rbindlist(all_spots)
spots[, species := factor(species, levels = c("Zebrafish", "Medaka"))]

cat(sprintf("\nTotal: %s spots loaded.\n\n", format(nrow(spots), big.mark = ",")))

# =============================================================================
# 1. NUCLEI COUNT PER FRAME
# =============================================================================

cat("── 1. Nuclei count per frame ──\n")

count_per_frame <- spots[, .(
  n_global   = .N,
  n_roi      = sum(IN_ROI),
  n_track_roi = sum(TRACK_IN_ROI)
), by = .(species, FRAME, time_min)]

# --- Normalise to frame-0 baseline for each species (fold change) ---
count_per_frame <- count_per_frame[order(species, FRAME)]
count_per_frame[, `:=`(
  n_global_t0   = n_global[1],
  n_roi_t0      = n_roi[1],
  n_track_roi_t0 = n_track_roi[1]
), by = species]
count_per_frame[, `:=`(
  fc_global    = n_global    / n_global_t0,
  fc_roi       = n_roi       / pmax(n_roi_t0, 1),
  fc_track_roi = n_track_roi / pmax(n_track_roi_t0, 1),
  frac_roi     = n_roi / n_global
)]

# --- 1a. Global count (faceted, free y -- each species on own scale) ---
p_count_global <- ggplot(count_per_frame, aes(x = time_min, y = n_global,
                                               color = species)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  scale_color_manual(values = species_colors) +
  scale_y_continuous(labels = label_comma()) +
  facet_wrap(~ species, scales = "free") +
  labs(title = "Total Nuclei Per Frame",
       subtitle = "Each species on its own scale",
       x = "Time (min)", y = "Number of nuclei") +
  theme_density() + theme(legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA))

# --- 1b. Global fold-change (shared axis, direct comparison) ---
p_count_fc <- ggplot(count_per_frame, aes(x = time_min, y = fc_global,
                                           color = species)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  scale_color_manual(values = species_colors) +
  labs(title = "Fold Change in Total Nuclei",
       subtitle = "Relative to frame 0 (shared axis)",
       x = "Time (min)", y = "Fold change (vs. t=0)") +
  theme_density()

# --- 1c. ROI count (faceted) ---
p_count_roi <- ggplot(count_per_frame, aes(x = time_min, y = n_roi,
                                            color = species)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  scale_color_manual(values = species_colors) +
  scale_y_continuous(labels = label_comma()) +
  facet_wrap(~ species, scales = "free") +
  labs(title = "Nuclei in Margin (ROI) Per Frame",
       subtitle = "Each species on its own scale",
       x = "Time (min)", y = "Number of nuclei in ROI") +
  theme_density() + theme(legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA))

# --- 1d. ROI fold-change ---
p_count_roi_fc <- ggplot(count_per_frame, aes(x = time_min, y = fc_roi,
                                              color = species)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  scale_color_manual(values = species_colors) +
  labs(title = "Fold Change in Margin Nuclei",
       subtitle = "Relative to frame 0 (shared axis)",
       x = "Time (min)", y = "Fold change (vs. t=0)") +
  theme_density()

# --- 1e. Fraction in ROI (already comparable) ---
p_frac_roi <- ggplot(count_per_frame, aes(x = time_min, y = frac_roi,
                                           color = species)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  scale_color_manual(values = species_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  facet_wrap(~ species, scales = "free_x") +
  labs(title = "Fraction of Nuclei in Margin",
       subtitle = "Proportion in ROI each frame",
       x = "Time (min)", y = "Fraction in ROI") +
  theme_density() + theme(
    strip.background = element_rect(fill = "grey95", color = NA))

p1 <- (p_count_global | p_count_fc) /
      (p_count_roi    | p_count_roi_fc) /
      p_frac_roi +
  plot_annotation(
    title    = "Nuclei Counts Per Frame -- Zebrafish vs. Medaka",
    subtitle = sprintf("ZF: %d frames (dt=%.0fs)  |  MK: %d frames (dt=%.0fs)",
                       uniqueN(all_spots$zebrafish$FRAME),
                       SPECIES$zebrafish$frame_interval_sec,
                       uniqueN(all_spots$medaka$FRAME),
                       SPECIES$medaka$frame_interval_sec),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey50", size = 10)
    )
  )

save_pdf(p1, "density_01_nuclei_count.pdf", width = 16, height = 14)

# =============================================================================
# 2. CUMULATIVE NUCLEI COUNT (UNIQUE SPOTS SEEN UP TO FRAME t)
# =============================================================================

cat("── 2. Cumulative nuclei ──\n")

# Cumulative: running total of the count per frame (integral of the curve)
# Since spots can appear and disappear across frames this is the running integral
# But more meaningfully: cumulative *unique* tracks or cumulative detections.
# Use cumulative sum of per-frame counts as area-under-curve measure.

cumul <- count_per_frame[order(species, FRAME)]
cumul[, cum_global := cumsum(n_global), by = species]
cumul[, cum_roi    := cumsum(n_roi),    by = species]

# Faceted absolute cumulative
p_cum_global <- ggplot(cumul, aes(x = time_min, y = cum_global, color = species)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = species_colors) +
  scale_y_continuous(labels = label_comma()) +
  facet_wrap(~ species, scales = "free") +
  labs(title = "Cumulative Nuclei Detections (Global)",
       subtitle = "Running sum -- each species on its own scale",
       x = "Time (min)", y = "Cumulative detections") +
  theme_density() + theme(legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA))

p_cum_roi <- ggplot(cumul, aes(x = time_min, y = cum_roi, color = species)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = species_colors) +
  scale_y_continuous(labels = label_comma()) +
  facet_wrap(~ species, scales = "free") +
  labs(title = "Cumulative Nuclei Detections (Margin)",
       subtitle = "Running sum -- each species on its own scale",
       x = "Time (min)", y = "Cumulative detections in ROI") +
  theme_density() + theme(legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA))

# Normalised cumulative (fold change from t=0)
cumul[, cum_global_fc := cum_global / cum_global[1], by = species]
cumul[, cum_roi_fc    := cum_roi    / pmax(cum_roi[1], 1), by = species]

p_cum_fc <- ggplot(cumul, aes(x = time_min, color = species)) +
  geom_line(aes(y = cum_global_fc), linewidth = 0.7) +
  geom_line(aes(y = cum_roi_fc), linewidth = 0.5, linetype = "longdash") +
  scale_color_manual(values = species_colors) +
  labs(title = "Cumulative Detections -- Fold Change",
       subtitle = "Solid = global, dashed = ROI  |  relative to t=0",
       x = "Time (min)", y = "Fold change (vs. t=0)") +
  theme_density()

p2 <- (p_cum_global | p_cum_roi) / p_cum_fc +
  plot_annotation(
    title = "Cumulative Nuclei Detections -- Zebrafish vs. Medaka",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
    )
  )

save_pdf(p2, "density_02_cumulative.pdf", width = 16, height = 10)

# =============================================================================
# 3. RATE OF CHANGE (Δn PER FRAME)
# =============================================================================

cat("── 3. Rate of change ──\n")

rate <- count_per_frame[order(species, FRAME)]
rate[, delta_global := n_global - shift(n_global, 1, type = "lag"), by = species]
rate[, delta_roi    := n_roi    - shift(n_roi,    1, type = "lag"), by = species]
rate <- rate[!is.na(delta_global)]

# Smooth with rolling mean (window = 5 frames)
smooth_window <- 5
rate[, delta_global_smooth := frollmean(delta_global, n = smooth_window,
                                         align = "center"), by = species]
rate[, delta_roi_smooth    := frollmean(delta_roi,    n = smooth_window,
                                         align = "center"), by = species]

# Faceted rate of change (free y per species)
p_rate_global <- ggplot(rate, aes(x = time_min)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_line(aes(y = delta_global, color = species), alpha = 0.15, linewidth = 0.3) +
  geom_line(aes(y = delta_global_smooth, color = species), linewidth = 0.7) +
  scale_color_manual(values = species_colors) +
  facet_wrap(~ species, scales = "free") +
  labs(title = "Rate of Change (Global)",
       subtitle = sprintf("delta-n per frame  |  smoothed: %d-frame rolling mean", smooth_window),
       x = "Time (min)", y = "Change in nuclei / frame") +
  theme_density() + theme(legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA))

p_rate_roi <- ggplot(rate, aes(x = time_min)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_line(aes(y = delta_roi, color = species), alpha = 0.15, linewidth = 0.3) +
  geom_line(aes(y = delta_roi_smooth, color = species), linewidth = 0.7) +
  scale_color_manual(values = species_colors) +
  facet_wrap(~ species, scales = "free") +
  labs(title = "Rate of Change (Margin / ROI)",
       subtitle = sprintf("delta-n per frame  |  smoothed: %d-frame rolling mean", smooth_window),
       x = "Time (min)", y = "Change in nuclei / frame in ROI") +
  theme_density() + theme(legend.position = "none",
    strip.background = element_rect(fill = "grey95", color = NA))

# Normalised rate: % change per frame relative to count at that frame
rate[, pct_change_global := delta_global / shift(n_global, 1, type = "lag") * 100,
     by = species]
rate[, pct_change_roi := delta_roi / shift(n_roi, 1, type = "lag") * 100,
     by = species]
rate[, pct_change_global_smooth := frollmean(pct_change_global, n = smooth_window,
                                              align = "center"), by = species]
rate[, pct_change_roi_smooth := frollmean(pct_change_roi, n = smooth_window,
                                           align = "center"), by = species]

p_rate_pct <- ggplot(rate, aes(x = time_min)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_line(aes(y = pct_change_global_smooth, color = species), linewidth = 0.7) +
  geom_line(aes(y = pct_change_roi_smooth, color = species),
            linewidth = 0.5, linetype = "longdash") +
  scale_color_manual(values = species_colors) +
  labs(title = "Relative Rate of Change (% per frame)",
       subtitle = "Solid = global, dashed = ROI  |  comparable across species",
       x = "Time (min)", y = "% change / frame") +
  theme_density()

p3 <- (p_rate_global | p_rate_roi) / p_rate_pct +
  plot_annotation(
    title = "Per-Frame Change in Nuclei Count -- Zebrafish vs. Medaka",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
    )
  )

save_pdf(p3, "density_03_rate_of_change.pdf", width = 16, height = 10)

# =============================================================================
# 4. GEODESIC SURFACE DENSITY
# =============================================================================
#
# Replicates embryo_viewer.py _compute_density():
#   1. Map (THETA_DEG, PHI_DEG) → unit-sphere XYZ
#   2. Adaptive radius: median k-th NN distance on 5000-point subsample
#      from the median-size frame
#   3. Per-frame: RANN::nn2 ball count within radius
#   4. Subtract self (−1)
#
# This gives per-spot neighbour count = local surface density proxy.
# We then aggregate: mean / median per frame → time series.
# =============================================================================

cat("── 4. Geodesic surface density ──\n")

compute_geodesic_density <- function(dt, k_nn = K_NN, subsample_n = SUBSAMPLE_N) {
  # Map to unit sphere
  theta_rad <- dt$THETA_DEG * pi / 180
  phi_rad   <- dt$PHI_DEG   * pi / 180

  xyz <- cbind(
    sin(theta_rad) * cos(phi_rad),
    sin(theta_rad) * sin(phi_rad),
    cos(theta_rad)
  )

  frames <- dt$FRAME
  unique_frames <- sort(unique(frames))

  # ── Adaptive radius from representative frame ──
  frame_sizes <- dt[, .N, by = FRAME]
  median_size <- median(frame_sizes$N)
  ref_frame   <- frame_sizes[which.min(abs(N - median_size)), FRAME]

  ref_idx  <- which(frames == ref_frame)
  ref_pts  <- xyz[ref_idx, , drop = FALSE]

  k_use <- min(k_nn, nrow(ref_pts) - 1)
  if (k_use < 2) {
    warning("Too few points for density; returning zeros")
    return(rep(0, nrow(dt)))
  }

  # Subsample
  set.seed(42)
  if (nrow(ref_pts) > subsample_n) {
    sub_idx <- sample.int(nrow(ref_pts), subsample_n)
    sub_pts <- ref_pts[sub_idx, , drop = FALSE]
  } else {
    sub_pts <- ref_pts
  }

  nn_res  <- nn2(sub_pts, sub_pts, k = k_use + 1)
  radius  <- median(nn_res$nn.dists[, k_use + 1])  # col k+1 = k-th NN (col 1 = self)
  radius  <- max(radius, 1e-6)

  cat(sprintf("    Adaptive radius = %.5f (unit-sphere chord)\n", radius))

  # ── Count neighbours per frame ──
  density <- numeric(nrow(dt))

  for (fr in unique_frames) {
    mask <- which(frames == fr)
    pts  <- xyz[mask, , drop = FALSE]
    n_fr <- nrow(pts)
    if (n_fr < 2) next

    # Use nn2 with searchtype = "radius" to count neighbours within radius.
    # RANN::nn2 doesn't support radius search directly.
    # Instead: query all pairs within radius using a k large enough, then count.
    # For efficiency, use the fact that max possible neighbours is bounded.
    # We'll use a reasonable k_max and count non-Inf distances.
    k_max <- min(n_fr - 1, max(100, as.integer(n_fr * 0.05)))

    nn <- nn2(pts, pts, k = min(k_max + 1, n_fr), searchtype = "radius",
              radius = radius)
    # nn$nn.dists: row i col j = distance to j-th NN. Self is col 1 (dist ≈ 0).
    # Cells beyond radius are Inf. Count non-Inf cols minus 1 (self).
    counts <- rowSums(nn$nn.dists <= radius) - 1L  # subtract self
    density[mask] <- counts
  }

  return(density)
}

# Compute density per species (heavy — ~minutes for medaka)
density_summaries <- list()

for (sp_key in names(SPECIES)) {
  sp <- SPECIES[[sp_key]]
  cat(sprintf("  Computing geodesic density for %s...\n", sp$name))

  dt <- all_spots[[sp_key]]

  # Global density
  cat("    [Global]\n")
  dt[, density_global := compute_geodesic_density(.SD)]

  # ROI-only density
  cat("    [Margin / ROI]\n")
  dt_roi <- dt[IN_ROI == TRUE]
  if (nrow(dt_roi) > 0) {
    dt_roi[, density_roi := compute_geodesic_density(.SD)]
    dt[IN_ROI == TRUE, density_roi := dt_roi$density_roi]
  }

  all_spots[[sp_key]] <- dt

  # Aggregate per frame
  agg_global <- dt[, .(
    mean_density   = mean(density_global, na.rm = TRUE),
    median_density = median(density_global, na.rm = TRUE),
    q25            = quantile(density_global, 0.25, na.rm = TRUE),
    q75            = quantile(density_global, 0.75, na.rm = TRUE),
    scope          = "Global"
  ), by = .(species, FRAME, time_min)]

  agg_roi <- dt[IN_ROI == TRUE, .(
    mean_density   = mean(density_roi, na.rm = TRUE),
    median_density = median(density_roi, na.rm = TRUE),
    q25            = quantile(density_roi, 0.25, na.rm = TRUE),
    q75            = quantile(density_roi, 0.75, na.rm = TRUE),
    scope          = "Margin (ROI)"
  ), by = .(species, FRAME, time_min)]

  density_summaries[[sp_key]] <- rbindlist(list(agg_global, agg_roi))

  cat(sprintf("    Done. Global density range: %.0f–%.0f neighbours\n",
              min(dt$density_global), max(dt$density_global)))
}

dens_all <- rbindlist(density_summaries)
dens_all[, species := factor(species, levels = c("Zebrafish", "Medaka"))]
dens_all[, scope   := factor(scope,   levels = c("Global", "Margin (ROI)"))]

# --- 4a. Mean geodesic density over time ---
p_dens_mean <- ggplot(dens_all, aes(x = time_min, color = species, fill = species)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
  geom_line(aes(y = median_density), linewidth = 0.7) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  facet_wrap(~ scope, scales = "free_y") +
  labs(title = "Geodesic Surface Density Over Time",
       subtitle = "Median neighbours within adaptive radius  |  ribbon = IQR",
       x = "Time (min)", y = "Neighbour count (geodesic)") +
  theme_density() +
  theme(strip.background = element_rect(fill = "grey95", color = NA))

save_pdf(p_dens_mean, "density_04_geodesic_density.pdf", width = 16, height = 7)

# --- 4b. Density distributions (violin) at selected time-points ---
# Pick 5 equally spaced frames per species
density_snapshots <- list()
for (sp_key in names(SPECIES)) {
  dt <- all_spots[[sp_key]]
  fr_range <- range(dt$FRAME)
  snap_frames <- round(seq(fr_range[1], fr_range[2], length.out = 5))

  snap <- dt[FRAME %in% snap_frames, .(species, FRAME, time_min,
                                         density_global, IN_ROI)]
  snap[, time_label := sprintf("%.0f min", time_min)]
  density_snapshots[[sp_key]] <- snap
}

snap_all <- rbindlist(density_snapshots)
snap_all[, species := factor(species, levels = c("Zebrafish", "Medaka"))]

p_dens_violin <- ggplot(snap_all, aes(x = reorder(time_label, FRAME),
                                       y = density_global, fill = species)) +
  geom_violin(alpha = 0.6, scale = "width", trim = TRUE, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8,
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = species_fills) +
  facet_wrap(~ species, scales = "free_x") +
  labs(title = "Density Distributions at Selected Time-Points",
       subtitle = "Geodesic neighbour count per nucleus",
       x = "Time-point", y = "Neighbours") +
  theme_density() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_pdf(p_dens_violin, "density_04b_density_violin.pdf", width = 14, height = 7)

# =============================================================================
# 5. COMBINED DASHBOARD
# =============================================================================

cat("── 5. Combined dashboard ──\n")

# Use normalised time (0–1) for direct comparison
for (sp_key in names(SPECIES)) {
  dt <- all_spots[[sp_key]]
  fr_range <- range(dt$FRAME)
  dt[, time_norm := (FRAME - fr_range[1]) / max(1, fr_range[2] - fr_range[1])]
  all_spots[[sp_key]] <- dt
}

# Normalised count
count_norm <- rbindlist(lapply(names(SPECIES), function(sp_key) {
  dt <- all_spots[[sp_key]]
  fr_range <- range(dt$FRAME)
  dt[, .(
    n_global     = .N,
    n_roi        = sum(IN_ROI),
    time_norm    = (FRAME[1] - fr_range[1]) / max(1, fr_range[2] - fr_range[1])
  ), by = .(species, FRAME)]
}))
count_norm[, species := factor(species, levels = c("Zebrafish", "Medaka"))]

# Normalised density
dens_norm <- dens_all[scope == "Global"]
dens_norm <- merge(dens_norm,
  rbindlist(lapply(names(SPECIES), function(sp_key) {
    dt <- all_spots[[sp_key]]
    fr_range <- range(dt$FRAME)
    unique(dt[, .(species, FRAME,
                  time_norm = (FRAME - fr_range[1]) / max(1, fr_range[2] - fr_range[1]))])
  })),
  by = c("species", "FRAME"), all.x = TRUE
)

# Merge fold-change into count_norm
count_norm <- count_norm[order(species, FRAME)]
count_norm[, `:=`(
  n_global_t0 = n_global[1],
  n_roi_t0    = n_roi[1]
), by = species]
count_norm[, `:=`(
  fc_global = n_global / n_global_t0,
  fc_roi    = n_roi    / pmax(n_roi_t0, 1)
)]

# Panel A: Fold-change total nuclei (shared axis, normalised time)
pA <- ggplot(count_norm, aes(x = time_norm, y = fc_global, color = species)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = species_colors) +
  labs(title = "A  Total Nuclei (fold change)",
       x = "Normalised time", y = "Fold change vs. t=0") +
  theme_density(base_size = 9) + theme(plot.title = element_text(size = 11))

# Panel B: Fold-change ROI nuclei
pB <- ggplot(count_norm, aes(x = time_norm, y = fc_roi, color = species)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = species_colors) +
  labs(title = "B  Margin Nuclei (fold change)",
       x = "Normalised time", y = "Fold change vs. t=0") +
  theme_density(base_size = 9) + theme(plot.title = element_text(size = 11))

# Panel C: Geodesic density (normalised time)
pC <- ggplot(dens_norm, aes(x = time_norm, color = species, fill = species)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.12, color = NA) +
  geom_line(aes(y = median_density), linewidth = 0.6) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  labs(title = "C  Geodesic Density", x = "Normalised time",
       y = "Median neighbours") +
  theme_density(base_size = 9) + theme(plot.title = element_text(size = 11))

# Panel D: Fraction in ROI (normalised time)
count_norm[, frac_roi := n_roi / n_global]
pD <- ggplot(count_norm, aes(x = time_norm, y = frac_roi, color = species)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = species_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "D  Fraction in Margin", x = "Normalised time",
       y = "% in ROI") +
  theme_density(base_size = 9) + theme(plot.title = element_text(size = 11))

p5 <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title    = "Density & Nuclei Count Dashboard -- Zebrafish vs. Medaka",
    subtitle = paste0(
      "ZF: ", format(nrow(all_spots$zebrafish), big.mark = ","), " spots, ",
      uniqueN(all_spots$zebrafish$FRAME), " frames, R = ",
      round(SPECIES$zebrafish$radius, 1), " um  |  ",
      "MK: ", format(nrow(all_spots$medaka), big.mark = ","), " spots, ",
      uniqueN(all_spots$medaka$FRAME), " frames, R = ",
      round(SPECIES$medaka$radius, 1), " um"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey50", size = 10)
    )
  )

save_pdf(p5, "density_05_combined_dashboard.pdf", width = 16, height = 10)

# =============================================================================
# 6. SEPARATE-PANEL PER SPECIES (side-by-side, shared y-axis where appropriate)
# =============================================================================

cat("── 6. Side-by-side species panels ──\n")

# Count per frame — faceted by species, free x (different durations)
p6_count <- ggplot(count_per_frame, aes(x = time_min)) +
  geom_area(aes(y = n_global), fill = "grey80", alpha = 0.5) +
  geom_area(aes(y = n_roi), fill = NA, color = "grey40", linewidth = 0.4,
            linetype = "dashed") +
  geom_line(aes(y = n_global, color = species), linewidth = 0.7) +
  geom_line(aes(y = n_roi, color = species), linewidth = 0.5, linetype = "longdash") +
  scale_color_manual(values = species_colors) +
  scale_y_continuous(labels = label_comma()) +
  facet_wrap(~ species, scales = "free") +
  labs(title = "Nuclei Count Over Time",
       subtitle = "Solid = global  |  Dashed = margin (ROI)",
       x = "Time (min)", y = "Nuclei") +
  theme_density() +
  theme(strip.background = element_rect(fill = "grey95", color = NA))

# Rate of change — faceted
p6_rate <- ggplot(rate, aes(x = time_min)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_line(aes(y = delta_global_smooth, color = species), linewidth = 0.7) +
  geom_line(aes(y = delta_roi_smooth, color = species), linewidth = 0.5,
            linetype = "longdash") +
  scale_color_manual(values = species_colors) +
  facet_wrap(~ species, scales = "free") +
  labs(title = "Rate of Change in Nuclei Count",
       subtitle = sprintf("Smoothed (%d-frame rolling mean)  |  Solid = global, Dashed = ROI",
                           smooth_window),
       x = "Time (min)", y = "Change in nuclei / frame") +
  theme_density() +
  theme(strip.background = element_rect(fill = "grey95", color = NA))

p6 <- p6_count / p6_rate +
  plot_annotation(
    title = "Per-Species Detail -- Zebrafish vs. Medaka",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
    )
  )

save_pdf(p6, "density_06_species_detail.pdf", width = 16, height = 12)

# =============================================================================
# DONE
# =============================================================================

cat("\n", strrep("═", 70), "\n")
cat("  All density plots saved to ", OUTPUT_DIR, "/\n")
cat(strrep("═", 70), "\n\n")
