# =============================================================================
# Zebrafish-Only Density & Flow Analysis
# =============================================================================
#
# Four publication-quality figures from the filtered+oriented zebrafish data:
#
#   1) NUCLEI COUNT  -- dual-axis: count + delta per frame
#   2) FLOW MAP      -- XY velocity arrows coloured by movement type (all frames)
#   3) DENSITY FOLD-CHANGE -- ROI vs outside geodesic density comparison
#   4) SHELL THINNING -- cross-section + thickness over time
#
# I/O directory: analysis_output_zebrafish_05112025/
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

IO_DIR  <- "analysis_output_zebrafish_05112025"

FRAME_INTERVAL_SEC <- 120
FRAME_INTERVAL_MIN <- FRAME_INTERVAL_SEC / 60

# Flow field (XY spatial bins) — poster-optimised, coherent with medaka
FLOW_BIN_XY  <- 30    # um
FLOW_MIN_N   <- 10    # per-direction minimum (split by movement type)
ARROW_SCALE  <- 60    # large arrows visible on poster
ARROW_HEAD   <- 0.18  # cm arrowhead
ARROW_LW     <- 1.2   # thick shaft for poster visibility

# Density vertical profile
Y_BIN_SIZE <- 15   # um

# Geodesic density (unit-sphere kNN)
K_NN        <- 16
SUBSAMPLE_N <- 5000

# =============================================================================
# HELPERS
# =============================================================================

theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      plot.subtitle    = element_text(size = base_size, hjust = 0.5, colour = "grey40"),
      strip.text       = element_text(face = "bold", size = base_size),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1)
    )
}

save_pdf <- function(plot, filename, width = 14, height = 8) {
  path <- file.path(IO_DIR, filename)
  ggsave(path, plot, width = width, height = height, device = "pdf")
  cat(sprintf("  -> Saved %s\n", path))
}

frame_to_min <- function(f) f * FRAME_INTERVAL_MIN

# =============================================================================
# DATA LOADING
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("  ZEBRAFISH DENSITY & FLOW ANALYSIS\n")
cat(strrep("=", 70), "\n\n")

cat("Loading oriented spots...\n")
spots <- fread(file.path(IO_DIR, "oriented_spots.csv"),
               select = c("TRACK_ID", "FRAME", "POSITION_T",
                           "POSITION_X", "POSITION_Y", "POSITION_Z",
                           "THETA_DEG", "PHI_DEG", "IN_ROI",
                           "RADIAL_VELOCITY", "RADIAL_VELOCITY_SMOOTH"))

MAX_FRAME <- spots[, max(FRAME)]

cat(sprintf("  Total: %s spots, frames %d-%d\n",
            format(nrow(spots), big.mark = ","),
            spots[, min(FRAME)], MAX_FRAME))

spots[, time_min := FRAME * FRAME_INTERVAL_MIN]

# =============================================================================
# 1. NUCLEI COUNT & DIVISION TIME (dual-axis single plot)
# =============================================================================

cat("\n-- 1. Nuclei count --\n")

cpf <- spots[, .(n_nuclei = .N), by = .(FRAME, time_min)]
setorder(cpf, FRAME)
cpf[, delta := n_nuclei - shift(n_nuclei, 1, type = "lag")]
cpf[, delta_sm := frollmean(delta, n = 11, align = "center")]

cat(sprintf("  Count range: %d - %d nuclei\n", min(cpf$n_nuclei), max(cpf$n_nuclei)))

# Net growth (percent)
growth_pct <- (cpf[FRAME == MAX_FRAME, n_nuclei] / cpf[FRAME == 0, n_nuclei] - 1) * 100
cat(sprintf("  Net growth: %.1f%% over %d frames\n", growth_pct, MAX_FRAME))

max_n     <- max(cpf$n_nuclei, na.rm = TRUE)
max_delta <- max(abs(cpf$delta_sm), na.rm = TRUE)
delta_sc  <- max_n / (max_delta * 2.5)

p1 <- ggplot(cpf, aes(x = time_min)) +
  geom_line(aes(y = n_nuclei), linewidth = 0.7, colour = "#B2182B") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70", linewidth = 0.3) +
  geom_line(aes(y = delta_sm * delta_sc), linewidth = 0.6, colour = "#2166AC", na.rm = TRUE) +
  scale_y_continuous(
    name   = "Nuclei count",
    labels = label_comma(),
    sec.axis = sec_axis(~ . / delta_sc, name = "Delta nuclei / frame (smoothed)")
  ) +
  annotate("text",
           x = max(cpf$time_min) * 0.55, y = max_n * 0.35,
           label = sprintf("Net growth: +%.1f%% over %d frames", growth_pct, MAX_FRAME),
           hjust = 0, size = 3.5, colour = "grey30") +
  labs(title = sprintf("Zebrafish -- Nuclei Count & Growth Rate (frames 0-%d)", MAX_FRAME),
       subtitle = "Red = nuclei count | Blue = per-frame change (11-frame smooth)",
       x = "Time (min)") +
  theme_pub() +
  theme(axis.title.y.left  = element_text(colour = "#B2182B"),
        axis.text.y.left   = element_text(colour = "#B2182B"),
        axis.title.y.right = element_text(colour = "#2166AC"),
        axis.text.y.right  = element_text(colour = "#2166AC"))

save_pdf(p1, "zf_density_01_nuclei_count.pdf", width = 14, height = 6)


# =============================================================================
# 2. FLOW VECTOR MAP  (XY space, Y reversed so AP on top -- all frames)
# =============================================================================

cat("\n-- 2. Flow vector map (all frames) --\n")

setorder(spots, TRACK_ID, FRAME)
spots[, `:=`(
  d_theta = THETA_DEG  - shift(THETA_DEG,  1, type = "lag"),
  dx      = POSITION_X - shift(POSITION_X, 1, type = "lag"),
  dy      = POSITION_Y - shift(POSITION_Y, 1, type = "lag"),
  dz      = POSITION_Z - shift(POSITION_Z, 1, type = "lag"),
  dt_f    = FRAME      - shift(FRAME,      1, type = "lag")
), by = TRACK_ID]

spots_v <- spots[!is.na(dx) & dt_f > 0]
spots_v[, `:=`(
  vx       = dx / dt_f,
  vy       = dy / dt_f,
  speed_3d = sqrt(dx^2 + dy^2 + dz^2) / dt_f
)]

# ---------------------------------------------------------------------------
# Per-point movement classification — mirrors embryo_viewer.py exactly
# (_compute_per_point_movement_type with half_window=20, drift_threshold=0.30)
#
# For each point in a track, look at net theta over a ±HALF_WIN window.
# Drift efficiency = |net_theta| / cumulative |step_theta|.
#   efficiency > 0.30 & net > 0  →  Downward (VP)   (red)
#   efficiency > 0.30 & net < 0  →  Upward   (AP)   (blue)
#   efficiency <= 0.30 or path < 0.5° → excluded
# ---------------------------------------------------------------------------
HALF_WIN        <- 20L
DRIFT_THRESHOLD <- 0.30

cat("  Computing per-point movement type (half_window=", HALF_WIN, ") ...\n")

setorder(spots, TRACK_ID, FRAME)
tid_vec   <- spots$TRACK_ID
theta_vec <- spots$THETA_DEG

n_all <- nrow(spots)
move_cat <- integer(n_all)  # 0=unclassified

boundaries <- which(diff(tid_vec) != 0)
starts <- c(1L, boundaries + 1L)
ends   <- c(boundaries, n_all)

for (k in seq_along(starts)) {
  s <- starts[k]; e <- ends[k]
  len <- e - s + 1L
  if (len < 3L) next

  t_seg  <- theta_vec[s:e]
  abs_dt <- abs(diff(t_seg))

  for (j in seq_len(len)) {
    lo <- max(1L, j - HALF_WIN)
    hi <- min(len, j + HALF_WIN + 1L)
    if ((hi - lo) < 2L) next

    net_theta  <- t_seg[hi] - t_seg[lo]
    path_theta <- sum(abs_dt[lo:(hi - 1L)])

    if (path_theta < 0.5) next

    eff <- abs(net_theta) / path_theta
    if (eff > DRIFT_THRESHOLD) {
      move_cat[s + j - 1L] <- if (net_theta > 0) 2L else 3L  # 2=down, 3=up
    }
  }
}
spots[, move_cat := move_cat]
spots[, move_label := fcase(
  move_cat == 2L, "Downward (VP)",
  move_cat == 3L, "Upward (AP)",
  default = NA_character_
)]

spots_v[, move_cat   := spots[!is.na(dx) & dt_f > 0, move_cat]]
spots_v[, move_label := spots[!is.na(dx) & dt_f > 0, move_label]]

cat(sprintf("  Classification (all frames): Downward %d, Upward %d, Unclassified %d\n",
            sum(spots_v$move_cat == 2L), sum(spots_v$move_cat == 3L),
            sum(spots_v$move_cat == 0L)))

# --- Two-panel flow map: Epiboly-only vs Epiboly + Ingression ---
# Ingression (upward movement) appears around frame 40.
EPOCH_EPIB <- c(0, 40)      # pure epiboly (only downward)
EPOCH_ING  <- c(40, 100)    # epiboly + ingression (upward/AP visible)

make_flow_map <- function(dt, epoch_name, epoch_frames,
                          show_types = c("Downward (VP)", "Upward (AP)")) {
  sub <- dt[FRAME >= epoch_frames[1] & FRAME < epoch_frames[2] &
            !is.na(move_label)]
  cat(sprintf("    %s: %s classified velocity pts\n", epoch_name,
              format(nrow(sub), big.mark = ",")))

  for (mt in c("Downward (VP)", "Upward (AP)")) {
    nn <- sum(sub$move_label == mt)
    cat(sprintf("      %s: %d (%.0f%%)\n", mt, nn, 100 * nn / nrow(sub)))
  }

  sub <- sub[move_label %in% show_types]

  fl_list <- list()
  for (mt in show_types) {
    s_mt <- sub[move_label == mt]
    if (nrow(s_mt) == 0L) next
    fl_mt <- s_mt[, .(
      mean_vx    = mean(vx, na.rm = TRUE),
      mean_vy    = mean(vy, na.rm = TRUE),
      mean_speed = mean(speed_3d, na.rm = TRUE),
      n = .N
    ), by = .(
      x_bin = floor(POSITION_X / FLOW_BIN_XY) * FLOW_BIN_XY + FLOW_BIN_XY / 2,
      y_bin = floor(POSITION_Y / FLOW_BIN_XY) * FLOW_BIN_XY + FLOW_BIN_XY / 2
    )]
    fl_mt <- fl_mt[n >= FLOW_MIN_N]
    fl_mt[, move_type := mt]
    fl_list <- c(fl_list, list(fl_mt))
  }
  fl <- rbindlist(fl_list, use.names = TRUE)
  fl[, move_type := factor(move_type,
        levels = c("Upward (AP)", "Downward (VP)"))]

  st <- fl[, .(n_bins = .N,
               pct = sprintf("%.0f%%", .N / nrow(fl) * 100)), by = move_type]
  setorder(st, move_type)
  cat(sprintf("      Bins: %s\n",
              paste0(st$move_type, ": ", st$n_bins, " (", st$pct, ")",
                     collapse = ", ")))

  cols <- c("Upward (AP)" = "#2166AC", "Downward (VP)" = "#B2182B")

  p <- ggplot(fl) +
    geom_segment(
      data = fl[move_type == "Downward (VP)"],
      aes(x = x_bin, y = y_bin,
          xend = x_bin + mean_vx * ARROW_SCALE,
          yend = y_bin + mean_vy * ARROW_SCALE,
          colour = move_type),
      arrow = arrow(length = unit(ARROW_HEAD, "cm"), type = "closed"),
      linewidth = ARROW_LW, alpha = 0.85) +
    geom_segment(
      data = fl[move_type == "Upward (AP)"],
      aes(x = x_bin, y = y_bin,
          xend = x_bin + mean_vx * ARROW_SCALE,
          yend = y_bin + mean_vy * ARROW_SCALE,
          colour = move_type),
      arrow = arrow(length = unit(ARROW_HEAD, "cm"), type = "closed"),
      linewidth = ARROW_LW, alpha = 0.95) +
    scale_colour_manual(values = cols, name = "Movement", drop = FALSE) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = epoch_name,
         subtitle = sprintf("Frames %d-%d (%.0f-%.0f min) | %d bins",
                            epoch_frames[1], epoch_frames[2],
                            frame_to_min(epoch_frames[1]),
                            frame_to_min(epoch_frames[2]),
                            nrow(fl)),
         x = expression("X (" * mu * "m)"),
         y = expression("Y (" * mu * "m)  [AP top, VP bottom]")) +
    theme_pub(base_size = 16) +
    theme(legend.text  = element_text(size = 15),
          legend.title = element_text(size = 16, face = "bold"),
          plot.title   = element_text(size = 18, face = "bold"),
          plot.subtitle = element_text(size = 13))
  p
}

p2a <- make_flow_map(spots_v, "Epiboly",
                     EPOCH_EPIB,
                     show_types = c("Downward (VP)"))
p2b <- make_flow_map(spots_v, "Epiboly + Ingression",
                     EPOCH_ING,
                     show_types = c("Downward (VP)"))

p2 <- (p2a | p2b) +
  plot_annotation(
    title = "Zebrafish -- XY Velocity Vector Field",
    subtitle = sprintf("Arrows = mean displacement \u00d7%d | bin = %d \u00b5m | Per-point rolling-window classification",
                       ARROW_SCALE, FLOW_BIN_XY),
    theme = theme(plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
                  plot.subtitle = element_text(size = 14, hjust = 0.5, colour = "grey40")))

save_pdf(p2, "zf_density_02_flow_maps.pdf", width = 22, height = 12)


# =============================================================================
# 3. ROI vs REST -- density redistribution (IN_ROI cells vs outside)
# =============================================================================

cat("\n-- 3. ROI vs Rest density comparison --\n")

WINDOW <- 5
roi_bounds <- fread(file.path(IO_DIR, "roi_bounds.csv"))
margin_y <- roi_bounds[parameter == "y_min", as.numeric(value)]

# ---------------------------------------------------------------------------
# Geodesic density (replicates embryo_viewer._compute_density)
# ---------------------------------------------------------------------------
compute_geodesic_density <- function(dt, k_nn = 16L, max_sub = 5000L) {
  theta_rad <- dt$THETA_DEG * pi / 180
  phi_rad   <- dt$PHI_DEG   * pi / 180
  coords <- cbind(
    sin(theta_rad) * cos(phi_rad),
    sin(theta_rad) * sin(phi_rad),
    cos(theta_rad)
  )
  frames <- dt$FRAME
  unique_frames <- sort(unique(frames))
  density <- numeric(nrow(dt))

  # Adaptive radius from a representative frame
  frame_sizes <- table(frames)
  median_size <- median(frame_sizes)
  ref_frame   <- as.integer(names(which.min(abs(frame_sizes - median_size))))
  ref_mask    <- frames == ref_frame
  ref_pts     <- coords[ref_mask, , drop = FALSE]

  k <- min(k_nn, nrow(ref_pts) - 1L)
  if (k < 2L) return(density)

  if (nrow(ref_pts) > max_sub) {
    set.seed(42)
    sub_idx <- sample(nrow(ref_pts), max_sub)
    sub_pts <- ref_pts[sub_idx, , drop = FALSE]
  } else {
    sub_pts <- ref_pts
  }

  nn_ref <- nn2(sub_pts, sub_pts, k = k + 1L)
  radius <- median(nn_ref$nn.dists[, k + 1L])
  radius <- max(radius, 1e-6)
  cat(sprintf("    Geodesic density: ref frame = %d, adaptive radius = %.5f\n",
              ref_frame, radius))

  for (fr in unique_frames) {
    mask    <- frames == fr
    pts_fr  <- coords[mask, , drop = FALSE]
    n_fr    <- sum(mask)
    if (n_fr < 2L) next
    k_search <- min(200L, n_fr)
    nn_fr    <- nn2(pts_fr, pts_fr, k = k_search)
    counts <- rowSums(nn_fr$nn.dists[, -1, drop = FALSE] <= radius)
    density[mask] <- counts
  }
  density
}

# Baseline (frame 0) and endpoint (frame 100) windows
ENDPOINT_FRAME <- 100
baseline <- spots[FRAME >= 0      & FRAME <= WINDOW]
endpoint <- spots[FRAME >= (ENDPOINT_FRAME - WINDOW) & FRAME <= ENDPOINT_FRAME]

n_roi_t0  <- nrow(baseline[IN_ROI == TRUE])
n_out_t0  <- nrow(baseline[IN_ROI == FALSE])
n_roi_end <- nrow(endpoint[IN_ROI == TRUE])
n_out_end <- nrow(endpoint[IN_ROI == FALSE])
cat(sprintf("  Frame 0:   IN_ROI = %d,  Outside = %d\n", n_roi_t0, n_out_t0))
cat(sprintf("  Frame %d: IN_ROI = %d,  Outside = %d\n", ENDPOINT_FRAME, n_roi_end, n_out_end))

# Compute geodesic density for baseline and endpoint
cat("  Computing geodesic density (baseline)...\n")
baseline[, geo_density := compute_geodesic_density(baseline)]
cat("  Computing geodesic density (endpoint)...\n")
endpoint[, geo_density := compute_geodesic_density(endpoint)]

mean_d_roi_t0  <- mean(baseline[IN_ROI == TRUE,  geo_density])
mean_d_out_t0  <- mean(baseline[IN_ROI == FALSE, geo_density])
mean_d_roi_end <- mean(endpoint[IN_ROI == TRUE,  geo_density])
mean_d_out_end <- mean(endpoint[IN_ROI == FALSE, geo_density])
cat(sprintf("  Mean density  ROI: %.1f -> %.1f (FC = %.2fx)\n",
            mean_d_roi_t0, mean_d_roi_end, mean_d_roi_end / mean_d_roi_t0))
cat(sprintf("  Mean density  Out: %.1f -> %.1f (FC = %.2fx)\n",
            mean_d_out_t0, mean_d_out_end, mean_d_out_end / mean_d_out_t0))

# --- Panel A: ZY cross-section, IN_ROI highlighted, t0 vs t_end ---
X_CENTRE <- median(spots$POSITION_X)
X_BAND   <- 30

set.seed(42)
n_sub <- 4000
cs_t0  <- baseline[abs(POSITION_X - X_CENTRE) <= X_BAND]
cs_end <- endpoint[abs(POSITION_X - X_CENTRE) <= X_BAND]
if (nrow(cs_t0)  > n_sub) cs_t0  <- cs_t0[sample(.N, n_sub)]
if (nrow(cs_end) > n_sub) cs_end <- cs_end[sample(.N, n_sub)]
cs_t0[,  epoch := "Frame 0"]
cs_end[, epoch := sprintf("Frame %d", ENDPOINT_FRAME)]
cs <- rbindlist(list(cs_t0, cs_end), fill = TRUE)
cs[, epoch := factor(epoch, levels = c("Frame 0", sprintf("Frame %d", ENDPOINT_FRAME)))]
cs[, region := fifelse(IN_ROI, "ROI", "Outside")]

p3a <- ggplot(cs, aes(x = POSITION_Z, y = POSITION_Y)) +
  geom_point(aes(colour = region), size = 0.6, alpha = 0.5) +
  scale_colour_manual(values = c("ROI" = "#D95F02", "Outside" = "grey55"),
                       name = "") +
  facet_wrap(~epoch) +
  scale_y_reverse() +
  coord_fixed() +
  labs(title = sprintf("ZY Cross-Section (X-slice \u00b1%d \u00b5m)", X_BAND),
       x = expression("Z (" * mu * "m)"),
       y = expression("Y (" * mu * "m)")) +
  theme_pub(base_size = 12) +
  theme(legend.position = "top")

# --- Panel B: Mean geodesic density (ROI vs Outside) ---
dens_dt <- data.table(
  region = rep(c("ROI", "Outside"), each = 2),
  epoch  = rep(c("Frame 0", sprintf("Frame %d", ENDPOINT_FRAME)), 2),
  density = c(mean_d_roi_t0, mean_d_roi_end, mean_d_out_t0, mean_d_out_end)
)
dens_dt[, epoch := factor(epoch, levels = c("Frame 0", sprintf("Frame %d", ENDPOINT_FRAME)))]

p3b <- ggplot(dens_dt, aes(x = region, y = density, fill = epoch)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f", density)),
            position = position_dodge(width = 0.6), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Frame 0" = "grey55",
                                setNames("#B2182B", sprintf("Frame %d", ENDPOINT_FRAME))),
                     name = "") +
  labs(title = "Mean Geodesic Density",
       subtitle = "Neighbours within adaptive radius (unit-sphere KNN)",
       x = NULL, y = "Mean neighbour count") +
  theme_pub(base_size = 12) +
  theme(legend.position = "top")

# --- Panel C: Density fold-change along Y, by region ---
MIN_BIN_COUNT <- 10

baseline[, y_bin := floor(POSITION_Y / Y_BIN_SIZE) * Y_BIN_SIZE + Y_BIN_SIZE / 2]
endpoint[, y_bin := floor(POSITION_Y / Y_BIN_SIZE) * Y_BIN_SIZE + Y_BIN_SIZE / 2]

# ROI cells
roi_t0  <- baseline[IN_ROI == TRUE,  .(mean_dens = mean(geo_density), n = .N), by = y_bin]
roi_end <- endpoint[IN_ROI == TRUE,  .(mean_dens = mean(geo_density), n = .N), by = y_bin]
fc_roi <- merge(roi_t0[n >= MIN_BIN_COUNT, .(y_bin, dens_t0  = mean_dens)],
                roi_end[n >= MIN_BIN_COUNT, .(y_bin, dens_end = mean_dens)],
                by = "y_bin", all = FALSE)
fc_roi[, log2_fc := log2(dens_end / pmax(dens_t0, 1e-6))]
fc_roi[, region := "ROI"]

# Outside cells
out_t0  <- baseline[IN_ROI == FALSE, .(mean_dens = mean(geo_density), n = .N), by = y_bin]
out_end <- endpoint[IN_ROI == FALSE, .(mean_dens = mean(geo_density), n = .N), by = y_bin]
fc_out <- merge(out_t0[n >= MIN_BIN_COUNT, .(y_bin, dens_t0  = mean_dens)],
                out_end[n >= MIN_BIN_COUNT, .(y_bin, dens_end = mean_dens)],
                by = "y_bin", all = FALSE)
fc_out[, log2_fc := log2(dens_end / pmax(dens_t0, 1e-6))]
fc_out[, region := "Outside"]

fc_all <- rbindlist(list(fc_roi, fc_out), fill = TRUE)
fc_all[log2_fc >  5, log2_fc :=  5]
fc_all[log2_fc < -5, log2_fc := -5]

p3c <- ggplot(fc_all, aes(x = y_bin, y = log2_fc, colour = region)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("ROI" = "#D95F02", "Outside" = "grey30"), name = "") +
  labs(title = sprintf("Density Fold Change (frame 0 vs %d)", ENDPOINT_FRAME),
       subtitle = "Mean geodesic density per Y-bin, log2 FC",
       x = expression("Y (" * mu * "m)"),
       y = "log2(fold change)") +
  theme_pub(base_size = 12) +
  theme(legend.position = "top")

p3 <- (p3a) / (p3b | p3c) +
  plot_layout(heights = c(1, 0.8)) +
  plot_annotation(
    title = sprintf("Zebrafish -- ROI vs Outside: Density Redistribution (frame 0 vs %d)", ENDPOINT_FRAME),
    subtitle = "Orange = ROI cells | Grey = outside ROI | Density = geodesic neighbour count (viewer method)",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey40")))

save_pdf(p3, "zf_density_03_density_fold_change.pdf", width = 16, height = 12)


# =============================================================================
# 4. SHELL THINNING -- Cross-section & thickness over time
# =============================================================================

cat("\n-- 4. Shell thinning analysis --\n")

spots_thick <- fread(file.path(IO_DIR, "oriented_spots.csv"),
                     select = c("TRACK_ID", "FRAME", "POSITION_X", "POSITION_Y",
                                "POSITION_Z", "RADIAL_DIST", "SPHERICAL_DEPTH",
                                "THETA_DEG", "PHI_DEG", "IN_ROI"))
spots_thick[, time_min := FRAME * FRAME_INTERVAL_MIN]

# Sphere params
sph <- fread(file.path(IO_DIR, "sphere_params.csv"))
SPH_R <- as.numeric(sph[parameter == "radius", value])
cat(sprintf("  Sphere radius: %.1f um\n", SPH_R))

# ---- Panel A: XZ cross-section at evenly spaced epochs --------------------
SPH_CX <- as.numeric(sph[parameter == "center_x", value])
SPH_CY <- as.numeric(sph[parameter == "center_y", value])
y_centre <- (margin_y + min(spots_thick$POSITION_Y)) / 2
Y_BAND   <- 15
X_BIN_THICK <- 15

# Spread epochs across the full timeline
epoch_frames <- round(seq(0, MAX_FRAME, length.out = 4))

cs_list <- lapply(epoch_frames, function(fr) {
  d <- spots_thick[FRAME >= (fr - WINDOW) & FRAME <= (fr + WINDOW) &
                     abs(POSITION_Y - y_centre) <= Y_BAND]
  d[, x_bin := floor(POSITION_X / X_BIN_THICK) * X_BIN_THICK + X_BIN_THICK / 2]
  d[, z_thickness := max(POSITION_Z) - min(POSITION_Z), by = x_bin]
  if (nrow(d) > 4000) d <- d[sample(.N, 4000)]
  d[, epoch := sprintf("Frame %d", fr)]
  d
})
cs <- rbindlist(cs_list, fill = TRUE)
cs[, epoch := factor(epoch, levels = sprintf("Frame %d", epoch_frames))]

thick_range <- quantile(cs$z_thickness, c(0.02, 0.98), na.rm = TRUE)

p4a <- ggplot(cs, aes(x = POSITION_X, y = POSITION_Z)) +
  geom_point(aes(colour = z_thickness), size = 0.5, alpha = 0.6) +
  scale_colour_viridis_c(option = "magma", name = "Z-thickness (\u00b5m)",
                          limits = thick_range, oob = scales::squish) +
  facet_wrap(~epoch, nrow = 1) +
  coord_fixed() +
  labs(title = "XZ Cross-Section -- Blastoderm Thickness",
       subtitle = sprintf("Thin Y-band (\u00b1%d \u00b5m) | colour = total Z-span per X bin",
                           Y_BAND),
       x = expression("X (" * mu * "m)"),
       y = expression("Z (" * mu * "m)")) +
  theme_pub(base_size = 11) +
  theme(legend.position = "bottom")

# ---- Panel B: X position vs fold change of Z-thickness --------------------
thick_baseline <- spots_thick[FRAME >= 0 & FRAME <= WINDOW &
                                abs(POSITION_Y - y_centre) <= Y_BAND]
thick_baseline[, x_bin := floor(POSITION_X / X_BIN_THICK) * X_BIN_THICK + X_BIN_THICK / 2]
zt0 <- thick_baseline[, .(z_thick_t0 = max(POSITION_Z) - min(POSITION_Z),
                           n = .N), by = x_bin]
zt0 <- zt0[n >= 5]

thick_endpoint <- spots_thick[FRAME >= (ENDPOINT_FRAME - WINDOW) & FRAME <= ENDPOINT_FRAME &
                                abs(POSITION_Y - y_centre) <= Y_BAND]
thick_endpoint[, x_bin := floor(POSITION_X / X_BIN_THICK) * X_BIN_THICK + X_BIN_THICK / 2]
zte <- thick_endpoint[, .(z_thick_end = max(POSITION_Z) - min(POSITION_Z),
                           n = .N), by = x_bin]
zte <- zte[n >= 5]

zfc <- merge(zt0[, .(x_bin, z_thick_t0)],
             zte[, .(x_bin, z_thick_end)],
             by = "x_bin", all = FALSE)
zfc[, log2_fc := log2(z_thick_end / pmax(z_thick_t0, 1e-3))]

cat(sprintf("  Z-thickness FC: %d X-bins, range log2FC = [%.2f, %.2f]\n",
            nrow(zfc), min(zfc$log2_fc), max(zfc$log2_fc)))

p4b <- ggplot(zfc, aes(x = x_bin, y = log2_fc)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.8, colour = "grey30") +
  geom_point(aes(fill = log2_fc), shape = 21, size = 3.5, colour = "grey20",
             stroke = 0.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 0, name = "log2(FC)") +
  labs(title = sprintf("Z-Thickness Fold Change (frame 0 vs %d)", ENDPOINT_FRAME),
       subtitle = "Blue = thinning, Red = thickening",
       x = expression("X (" * mu * "m)"),
       y = "log2(fold change)") +
  theme_pub(base_size = 11)

# ---- Compose Figure 4 ---------------------------------------------------
p4 <- p4a / p4b +
  plot_layout(heights = c(1, 0.8)) +
  plot_annotation(
    title = "Zebrafish -- Blastoderm Thinning During Epiboly",
    subtitle = sprintf("Z-span per X column | Y-band \u00b1%d \u00b5m | %d \u00b5m X bins",
                        Y_BAND, X_BIN_THICK),
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5, colour = "grey40")))

save_pdf(p4, "zf_density_04_shell_thinning.pdf", width = 16, height = 10)


# =============================================================================
# DONE
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("  ZEBRAFISH DENSITY ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n")
cat(sprintf("  Output: %s/zf_density_*.pdf\n\n", IO_DIR))
