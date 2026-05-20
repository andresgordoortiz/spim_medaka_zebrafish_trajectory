# =============================================================================
# ONE-SHOT MEDAKA vs ZEBRAFISH COMPARISON
# =============================================================================
# Generates the plot set requested:
#   (1) margin cell density over time + fold difference (medaka/zebrafish)
#   (2) movement direction over time (epiboly fraction, mean v_epiboly,
#       per-step type composition)
#   (3) thickness in same plot (both species)
#   (4) flow fields per hour, both species, + majority movement type per hour
#   (5) matched-track comparison (>=20 µm displacement, 10-60 min duration):
#       straightness, mean turning angle, MSD
#
# Uses the same oriented-track inputs as gastrulation_dynamics_comparison.R
# but writes to a separate output folder so existing PDFs are untouched.
# =============================================================================

suppressPackageStartupMessages({
  source("renv/activate.R")
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(viridis)
  library(readr)
})

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

OUT_DIR <- "analysis_output_oneshot"
dir.create(OUT_DIR, showWarnings = FALSE)

MEDAKA_INPUT  <- "oriented_medaka_ultrack"
ZEBRAFISH_INPUT <- "oriented_zebrafish_ultrack"

# Output dirs holding sphere_params, landmarks and thickness CSVs
MEDAKA_OUT_LEGACY <- "analysis_output_medaka"        # has medaka_thickness_data.csv
MEDAKA_OUT_NEW    <- "analysis_output_medaka_25082025"
ZEB_OUT_LEGACY    <- "analysis_output_zebrafish"
ZEB_OUT_NEW       <- "analysis_output_zebrafish_05112025"

MEDAKA_FI  <- 30   # frame interval (s)
ZEB_FI     <- 120

MEDAKA_VOXEL_UM <- 1.05152
ZEB_VOXEL_UM    <- 1.24785

MEDAKA_MIN_FRAMES <- 20L   # 10 min
ZEB_MIN_FRAMES    <- 5L    # 10 min

MEDAKA_VEL_LAG    <- 4L    # 2 min
ZEB_VEL_LAG       <- 1L
MEDAKA_SMOOTH_K   <- 5L
ZEB_SMOOTH_K      <- 3L

MARGIN_WIDTH_DEG  <- 5     # band ±5° around margin landmark

FLOW_BIN_UM       <- 30
FLOW_MIN_N        <- 10
ARROW_SCALE       <- 60
ARROW_HEAD        <- 0.14
ARROW_LW          <- 0.7
VORT_SWIRL_THRESH <- 0.3

species_colors <- c("Medaka" = "#D73027", "Zebrafish" = "#4575B4")
step_type_cols <- c("Epiboly" = "#B2182B", "Animalward" = "#2166AC",
                    "Convergence" = "#7B3294")
flow_cols <- c("Upward (AP)" = "#2166AC", "Circular" = "#1A9850",
               "Downward (VP)" = "#B2182B")

theme_pub <- function(bs = 11) {
  theme_minimal(base_size = bs) +
    theme(
      plot.title    = element_text(face = "bold", size = bs + 1, hjust = 0),
      plot.subtitle = element_text(size = bs - 1, color = "grey40", hjust = 0),
      strip.text    = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey96", color = NA),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
}

save_pdf <- function(p, name, w = 12, h = 8) {
  path <- file.path(OUT_DIR, name)
  ok <- tryCatch({
    ggsave(path, p, width = w, height = h, device = cairo_pdf)
    TRUE
  }, error = function(e) {
    message("  cairo_pdf failed (", conditionMessage(e), "); falling back to pdf")
    ggsave(path, p, width = w, height = h, device = "pdf")
    FALSE
  })
  if (!file.exists(path)) {
    message("  WARNING: file not written, retrying with base pdf device")
    ggsave(path, p, width = w, height = h, device = "pdf")
  }
  cat(sprintf("  saved %s\n", path))
}

banner <- function(x) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("  ", x, "\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")
}

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------

banner("LOAD INPUT DATA")

cat("Loading raw oriented tracks...\n")
sp_m <- fread(file.path(MEDAKA_INPUT, "oriented_tracks_medaka.csv"),
              showProgress = FALSE)
sp_z <- fread(file.path(ZEBRAFISH_INPUT, "oriented_tracks_zebrafish.csv"),
              showProgress = FALSE)

# Voxel calibration → microns
for (col in c("POSITION_X", "POSITION_Y", "POSITION_Z", "RADIAL_DIST",
              "SPHERICAL_DEPTH")) {
  sp_m[[col]] <- sp_m[[col]] * MEDAKA_VOXEL_UM
  sp_z[[col]] <- sp_z[[col]] * ZEB_VOXEL_UM
}

# Filter short tracks
filter_short <- function(dt, min_n) {
  ids <- dt[, .N, by = TRACK_ID][N >= min_n, TRACK_ID]
  dt[TRACK_ID %in% ids]
}
sp_m <- filter_short(sp_m, MEDAKA_MIN_FRAMES)
sp_z <- filter_short(sp_z, ZEB_MIN_FRAMES)

sp_m[, `:=`(species = "Medaka",    time_min = FRAME * MEDAKA_FI / 60)]
sp_z[, `:=`(species = "Zebrafish", time_min = FRAME * ZEB_FI    / 60)]

cat(sprintf("  Medaka:    %s spots / %s tracks (frames %d-%d, %.0f min)\n",
            format(nrow(sp_m), big.mark = ","),
            format(uniqueN(sp_m$TRACK_ID), big.mark = ","),
            min(sp_m$FRAME), max(sp_m$FRAME), max(sp_m$time_min)))
cat(sprintf("  Zebrafish: %s spots / %s tracks (frames %d-%d, %.0f min)\n",
            format(nrow(sp_z), big.mark = ","),
            format(uniqueN(sp_z$TRACK_ID), big.mark = ","),
            min(sp_z$FRAME), max(sp_z$FRAME), max(sp_z$time_min)))

# Sphere params and landmarks
read_sphere <- function(dir, vox) {
  d <- fread(file.path(dir, "sphere_params.csv"))
  as.numeric(d[parameter == "radius", value]) * vox
}
R_M <- read_sphere(MEDAKA_INPUT, MEDAKA_VOXEL_UM)
R_Z <- read_sphere(ZEBRAFISH_INPUT, ZEB_VOXEL_UM)
cat(sprintf("  Sphere radius: medaka %.1f µm, zebrafish %.1f µm\n", R_M, R_Z))

read_margin <- function(dir) {
  d <- fread(file.path(dir, "gastrulation_landmarks.csv"))
  as.numeric(d[landmark == "margin", theta])
}
MARGIN_M <- read_margin(MEDAKA_INPUT)
MARGIN_Z <- read_margin(ZEBRAFISH_INPUT)
cat(sprintf("  Margin θ: medaka %.2f°, zebrafish %.2f°\n", MARGIN_M, MARGIN_Z))

# Estimate covered phi range (azimuth coverage of the imaged half) per species
phi_range <- function(dt) {
  q <- quantile(dt$PHI_DEG, c(0.02, 0.98), na.rm = TRUE)
  as.numeric(diff(q)) * pi / 180
}
PHI_RNG_M <- phi_range(sp_m)
PHI_RNG_Z <- phi_range(sp_z)
cat(sprintf("  Phi coverage: medaka %.2f rad, zebrafish %.2f rad\n",
            PHI_RNG_M, PHI_RNG_Z))

# -----------------------------------------------------------------------------
# Velocity / per-step metrics (smoothed)
# -----------------------------------------------------------------------------

banner("COMPUTE PER-STEP VELOCITIES")

compute_vel <- function(df, fi_sec, vel_lag, smooth_k) {
  fi_min <- fi_sec / 60
  dt <- copy(df); setkey(dt, TRACK_ID, FRAME)
  for (col in c("POSITION_X", "POSITION_Y", "POSITION_Z",
                "RADIAL_DIST", "THETA_DEG", "PHI_DEG")) {
    sm <- paste0(col, "_SM")
    dt[, (sm) := frollmean(get(col), n = smooth_k, align = "center"),
       by = TRACK_ID]
    dt[is.na(get(sm)), (sm) := get(col)]
  }
  dt[, `:=`(
    dx = POSITION_X_SM - shift(POSITION_X_SM, vel_lag),
    dy = POSITION_Y_SM - shift(POSITION_Y_SM, vel_lag),
    dz = POSITION_Z_SM - shift(POSITION_Z_SM, vel_lag),
    dt_f = FRAME - shift(FRAME, vel_lag),
    dTheta = THETA_DEG_SM - shift(THETA_DEG_SM, vel_lag),
    dPhi   = PHI_DEG_SM   - shift(PHI_DEG_SM,   vel_lag),
    lag_dx_prev = shift(POSITION_X_SM, vel_lag) -
      shift(POSITION_X_SM, 2 * vel_lag),
    lag_dy_prev = shift(POSITION_Y_SM, vel_lag) -
      shift(POSITION_Y_SM, 2 * vel_lag),
    lag_dz_prev = shift(POSITION_Z_SM, vel_lag) -
      shift(POSITION_Z_SM, 2 * vel_lag)
  ), by = TRACK_ID]
  dt[, dPhi := fifelse(abs(dPhi) > 180, dPhi - sign(dPhi) * 360, dPhi)]
  dt[, `:=`(
    disp_3d   = sqrt(dx^2 + dy^2 + dz^2),
    theta_rad = THETA_DEG * pi / 180
  )]
  dt[, `:=`(
    inst_speed    = disp_3d / (dt_f * fi_min),
    vx_um_min     = dx / (dt_f * fi_min),
    vy_um_min     = dy / (dt_f * fi_min),
    v_epiboly     = RADIAL_DIST * (dTheta * pi / 180) / (dt_f * fi_min),
    v_convergence = RADIAL_DIST * sin(theta_rad) *
                    (dPhi * pi / 180) / (dt_f * fi_min)
  )]
  # turning angle vs previous step
  dt[, `:=`(
    mag_prev = sqrt(lag_dx_prev^2 + lag_dy_prev^2 + lag_dz_prev^2),
    dot_prev = dx * lag_dx_prev + dy * lag_dy_prev + dz * lag_dz_prev
  )]
  dt[, cos_a := dot_prev / (disp_3d * mag_prev)]
  dt[, turning_angle := acos(pmin(pmax(cos_a, -1), 1)) * 180 / pi]
  dt[, c("lag_dx_prev","lag_dy_prev","lag_dz_prev","mag_prev","dot_prev",
         "cos_a") := NULL]
  dt[!is.na(dt_f) & dt_f == vel_lag]
}

vel_m <- compute_vel(sp_m, MEDAKA_FI, MEDAKA_VEL_LAG, MEDAKA_SMOOTH_K)
vel_z <- compute_vel(sp_z, ZEB_FI,    ZEB_VEL_LAG,    ZEB_SMOOTH_K)

vel_m[, step_type := fcase(
  abs(v_convergence) > abs(v_epiboly), "Convergence",
  v_epiboly > 0, "Epiboly",
  default = "Animalward"
)]
vel_z[, step_type := fcase(
  abs(v_convergence) > abs(v_epiboly), "Convergence",
  v_epiboly > 0, "Epiboly",
  default = "Animalward"
)]
vel_m[, species := "Medaka"];  vel_z[, species := "Zebrafish"]
vel_m[, time_min := FRAME * MEDAKA_FI / 60]
vel_z[, time_min := FRAME * ZEB_FI    / 60]

cat(sprintf("  velocity steps: medaka %s, zebrafish %s\n",
            format(nrow(vel_m), big.mark = ","),
            format(nrow(vel_z), big.mark = ",")))

# =============================================================================
# FIGURE 1: MARGIN CELL DENSITY OVER TIME + FOLD DIFFERENCE
# =============================================================================

banner("FIGURE 1 — margin cell density over time")

# Surface density at the margin band [margin - W, margin + W]
# area_band = R² * Δθ * sin(θ_mid) * Δφ  (Δφ from observed phi coverage)
margin_density <- function(dt, R, margin_theta, phi_rng, fi_sec,
                            width = MARGIN_WIDTH_DEG) {
  lo <- margin_theta - width; hi <- margin_theta + width
  band <- dt[THETA_DEG >= lo & THETA_DEG <= hi]
  if (nrow(band) == 0) return(data.table())
  dth_rad  <- (hi - lo) * pi / 180
  thm_rad  <- margin_theta * pi / 180
  area_um2 <- R^2 * dth_rad * sin(thm_rad) * phi_rng
  band[, .(n = .N, time_min = FRAME[1] * fi_sec / 60), by = FRAME][,
      .(FRAME, time_min, n_margin = n,
        density_per_1000um2 = n / area_um2 * 1000,
        area_um2 = area_um2)]
}

md_m <- margin_density(sp_m, R_M, MARGIN_M, PHI_RNG_M, MEDAKA_FI)
md_z <- margin_density(sp_z, R_Z, MARGIN_Z, PHI_RNG_Z, ZEB_FI)
md_m[, species := "Medaka"];  md_z[, species := "Zebrafish"]

# Bin by 30 min and compute fold difference medaka/zebrafish
BIN_MIN <- 30
md_m[, time_bin := floor(time_min / BIN_MIN) * BIN_MIN + BIN_MIN / 2]
md_z[, time_bin := floor(time_min / BIN_MIN) * BIN_MIN + BIN_MIN / 2]
binned <- rbind(md_m, md_z)[, .(
  density = mean(density_per_1000um2, na.rm = TRUE),
  n_frames = .N
), by = .(species, time_bin)]

fold <- dcast(binned, time_bin ~ species, value.var = "density")
setnames(fold, c("Medaka", "Zebrafish"), c("med", "zeb"), skip_absent = TRUE)
fold[, fold_med_over_zeb := med / zeb]

cat("  Margin band density summary (per 1000 µm²):\n")
print(rbind(md_m, md_z)[, .(min = min(density_per_1000um2),
                              median = median(density_per_1000um2),
                              max = max(density_per_1000um2)),
                          by = species])
cat(sprintf("  Mean fold (Medaka / Zebrafish) over time: %.2f×\n",
            mean(fold$fold_med_over_zeb, na.rm = TRUE)))

# Save tidy CSV
fwrite(rbind(md_m, md_z), file.path(OUT_DIR, "margin_density_per_frame.csv"))
fwrite(binned, file.path(OUT_DIR, "margin_density_binned.csv"))
fwrite(fold, file.path(OUT_DIR, "margin_density_fold.csv"))

p1a <- ggplot(rbind(md_m, md_z),
              aes(time_min, density_per_1000um2, color = species)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_smooth(se = TRUE, method = "loess", span = 0.3, linewidth = 1) +
  scale_color_manual(values = species_colors) +
  labs(title = "Cell density in margin band (±5° θ) over time",
       subtitle = sprintf("Margin θ: medaka %.1f°, zebrafish %.1f° — surface density",
                          MARGIN_M, MARGIN_Z),
       x = "Time (min)", y = expression("nuclei / 1000 " * mu * "m" ^ 2),
       color = NULL) +
  theme_pub()

p1b <- ggplot(binned, aes(time_bin, density, color = species)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_color_manual(values = species_colors) +
  labs(title = sprintf("Mean margin density (binned %d-min)", BIN_MIN),
       x = "Time (min)", y = expression("nuclei / 1000 " * mu * "m" ^ 2),
       color = NULL) +
  theme_pub()

p1c <- ggplot(fold[!is.na(fold_med_over_zeb)],
              aes(time_bin, fold_med_over_zeb)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1, color = "#444444") +
  geom_point(size = 2, color = "#444444") +
  geom_text(aes(label = sprintf("%.2f×", fold_med_over_zeb)),
            vjust = -0.6, size = 3) +
  labs(title = "Fold difference (Medaka / Zebrafish)",
       subtitle = "Ratio of mean margin densities in 30-min bins",
       x = "Time (min)", y = "Medaka density ÷ Zebrafish density") +
  theme_pub()

fig1 <- (p1a / p1b / p1c) +
  plot_annotation(
    title = "Margin cell density over time — Medaka vs Zebrafish",
    theme = theme(plot.title = element_text(face = "bold", size = 14)))
save_pdf(fig1, "01_margin_density.pdf", w = 12, h = 14)

# =============================================================================
# FIGURE 2: MOVEMENT DIRECTION OVER TIME (both species)
# =============================================================================

banner("FIGURE 2 — movement direction over time")

TIME_BIN <- 30
vel_all <- rbind(
  vel_m[, .(species, time_min, v_epiboly, v_convergence, step_type)],
  vel_z[, .(species, time_min, v_epiboly, v_convergence, step_type)]
)
vel_all[, time_bin := floor(time_min / TIME_BIN) * TIME_BIN + TIME_BIN / 2]

dir_summary <- vel_all[!is.na(v_epiboly),
  .(epi_pct        = 100 * mean(v_epiboly > 0),
    mean_v_epi     = mean(v_epiboly),
    mean_v_conv    = mean(v_convergence, na.rm = TRUE),
    abs_v_conv     = mean(abs(v_convergence), na.rm = TRUE),
    n = .N),
  by = .(species, time_bin)][n >= 50]

p2a <- ggplot(dir_summary, aes(time_bin, epi_pct, color = species)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_color_manual(values = species_colors) +
  labs(title = "Fraction of vegetalward (epiboly) steps over time",
       x = "Time (min)", y = "% steps with v_epiboly > 0", color = NULL) +
  theme_pub()

p2b <- ggplot(dir_summary, aes(time_bin, mean_v_epi, color = species)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_color_manual(values = species_colors) +
  labs(title = "Mean epiboly velocity over time",
       x = "Time (min)", y = expression("v"["epiboly"] * " (" * mu * "m/min)"),
       color = NULL) +
  theme_pub()

p2c <- ggplot(dir_summary, aes(time_bin, abs_v_conv, color = species)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_color_manual(values = species_colors) +
  labs(title = "|Convergence velocity| over time",
       x = "Time (min)",
       y = expression("|v"["conv"] * "| (" * mu * "m/min)"), color = NULL) +
  theme_pub()

step_comp <- vel_all[!is.na(step_type),
  .(n = .N), by = .(species, time_bin, step_type)]
step_comp[, pct := 100 * n / sum(n), by = .(species, time_bin)]
step_comp[, step_type := factor(step_type,
                                 levels = c("Epiboly", "Animalward",
                                            "Convergence"))]

p2d <- ggplot(step_comp, aes(time_bin, pct, fill = step_type)) +
  geom_area(alpha = 0.9, color = "white", linewidth = 0.2) +
  facet_wrap(~ species, ncol = 1) +
  scale_fill_manual(values = step_type_cols) +
  labs(title = "Per-step direction composition (stacked %)",
       x = "Time (min)", y = "% of steps", fill = NULL) +
  theme_pub()

fig2 <- (p2a + p2b) / (p2c + p2d) +
  plot_annotation(title = "Movement direction over time",
                  theme = theme(plot.title = element_text(face = "bold",
                                                          size = 14)))
save_pdf(fig2, "02_movement_direction.pdf", w = 14, h = 11)

fwrite(dir_summary, file.path(OUT_DIR, "direction_summary.csv"))
fwrite(step_comp, file.path(OUT_DIR, "step_type_composition.csv"))

# =============================================================================
# FIGURE 3: THICKNESS — both species in same plot
# Top:    thickness metric (IQR of SPHERICAL_DEPTH) per zone over time
# Bottom: depth distribution that the IQR summarises
#         (p10–p90 outer ribbon, p25–p75 inner ribbon = the IQR, median line)
# =============================================================================

banner("FIGURE 3 — thickness")

# Zone widths (match per-species pipelines)
ZONE_MARG_W_M  <- 5;  ZONE_PREM_W_M <- 10   # medaka
ZONE_MARG_W_Z  <- 8;  ZONE_PREM_W_Z <- 15   # zebrafish

assign_zone <- function(theta, margin_theta, marg_w, prem_w) {
  zm_top <- margin_theta - marg_w
  zp_top <- zm_top       - prem_w
  factor(fcase(
    theta > margin_theta, "Sub-marginal",
    theta > zm_top,       "Margin layer",
    theta > zp_top,       "Pre-marginal",
    default               = "Animal cap"),
    levels = c("Animal cap", "Pre-marginal", "Margin layer", "Sub-marginal"))
}

sp_m[, zone := assign_zone(THETA_DEG, MARGIN_M, ZONE_MARG_W_M, ZONE_PREM_W_M)]
sp_z[, zone := assign_zone(THETA_DEG, MARGIN_Z, ZONE_MARG_W_Z, ZONE_PREM_W_Z)]

# Per-cell depth distribution in 10-min bins per zone, per species
depth_dist <- rbind(
  sp_m[!is.na(zone) & is.finite(SPHERICAL_DEPTH),
       .(time_bin = floor(time_min / 10) * 10,
         zone, depth_um = SPHERICAL_DEPTH, species = "Medaka")],
  sp_z[!is.na(zone) & is.finite(SPHERICAL_DEPTH),
       .(time_bin = floor(time_min / 10) * 10,
         zone, depth_um = SPHERICAL_DEPTH, species = "Zebrafish")]
)
depth_dist[, species := factor(species, levels = c("Medaka", "Zebrafish"))]

depth_summary <- depth_dist[,
  .(p10 = quantile(depth_um, 0.10),
    p25 = quantile(depth_um, 0.25),
    p50 = quantile(depth_um, 0.50),
    p75 = quantile(depth_um, 0.75),
    p90 = quantile(depth_um, 0.90),
    thickness_iqr = IQR(depth_um),
    n_cells = .N),
  by = .(species, zone, time_bin)][n_cells >= 20]

# Top panel — the metric used to track thickness
p3a <- ggplot(depth_summary, aes(time_bin, thickness_iqr, color = species)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.4, alpha = 0.7) +
  facet_wrap(~ zone, nrow = 1, scales = "free_x") +
  scale_color_manual(values = species_colors) +
  labs(title = "Tissue thickness — IQR of cell depth per zone",
       x = "Time (min)", y = "thickness IQR (µm)", color = NULL) +
  theme_pub()

# Bottom panel — the underlying depth distribution that defines the IQR
p3b <- ggplot(depth_summary, aes(x = time_bin)) +
  geom_ribbon(aes(ymin = p10, ymax = p90, fill = species), alpha = 0.18) +
  geom_ribbon(aes(ymin = p25, ymax = p75, fill = species), alpha = 0.40) +
  geom_line(aes(y = p50, color = species), linewidth = 0.9) +
  facet_grid(species ~ zone, scales = "free") +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_color_manual(values = species_colors, guide = "none") +
  scale_y_reverse() +    # depth: 0 = surface, positive = deeper
  labs(title = "Depth distribution underlying the thickness measurement",
       subtitle = "outer band = p10–p90, inner band = IQR (p25–p75), line = median; y-axis reversed (surface up)",
       x = "Time (min)", y = "spherical depth (µm, surface up)") +
  theme_pub()

fig3 <- (p3a / p3b) +
  plot_layout(heights = c(1, 1.8)) +
  plot_annotation(
    title = "Tissue thickness over time — Medaka vs Zebrafish",
    subtitle = "Top = the IQR metric; bottom = the cell-depth cloud it summarises",
    theme = theme(plot.title = element_text(face = "bold", size = 14)))
save_pdf(fig3, "03_thickness.pdf", w = 16, h = 12)

fwrite(depth_summary, file.path(OUT_DIR, "thickness_depth_distribution.csv"))

# Quick numeric check — does the margin layer actually thicken?
margin_thick <- depth_summary[zone == "Margin layer"]
for (sp in c("Medaka", "Zebrafish")) {
  dd <- margin_thick[species == sp][order(time_bin)]
  if (nrow(dd) >= 3) {
    t0 <- dd$thickness_iqr[1]
    tmax <- max(dd$thickness_iqr)
    tend <- tail(dd$thickness_iqr, 1)
    tmax_t <- dd$time_bin[which.max(dd$thickness_iqr)]
    cat(sprintf("  %-9s margin IQR: t=0 %.1f µm | peak %.1f µm at t=%g min | end %.1f µm (Δpeak=%+.1f%%, Δend=%+.1f%%)\n",
                sp, t0, tmax, tmax_t, tend,
                100 * (tmax - t0) / t0, 100 * (tend - t0) / t0))
  }
}

# =============================================================================
# FIGURE 4: FLOW FIELDS PER HOUR + MAJORITY MOVEMENT TYPE
# =============================================================================

banner("FIGURE 4 — hourly flow fields")

# Add hourly bin
vel_m[, hour := floor(time_min / 60)]
vel_z[, hour := floor(time_min / 60)]

build_hourly_flow <- function(vel, sp_label, bin_um = FLOW_BIN_UM,
                              min_n = FLOW_MIN_N,
                              swirl_thr = VORT_SWIRL_THRESH) {
  bs <- bin_um
  res <- list()
  for (h in sort(unique(vel$hour))) {
    sub <- vel[hour == h]
    if (nrow(sub) < 200) next
    fl <- sub[, .(
      mean_vx = mean(vx_um_min, na.rm = TRUE),
      mean_vy = mean(vy_um_min, na.rm = TRUE),
      mean_speed = mean(inst_speed, na.rm = TRUE),
      n = .N
    ), by = .(
      x_bin = floor(POSITION_X / bs) * bs + bs / 2,
      y_bin = floor(POSITION_Y / bs) * bs + bs / 2
    )][n >= min_n]
    if (nrow(fl) < 5) next

    # 3x3 smoothing
    setkey(fl, x_bin, y_bin)
    svx <- numeric(nrow(fl)); svy <- numeric(nrow(fl))
    for (i in seq_len(nrow(fl))) {
      xb <- fl$x_bin[i]; yb <- fl$y_bin[i]
      vxs <- c(); vys <- c()
      for (dx_ in c(-bs, 0, bs)) for (dy_ in c(-bs, 0, bs)) {
        nb <- fl[.(xb + dx_, yb + dy_)]
        if (nrow(nb) == 1L) { vxs <- c(vxs, nb$mean_vx); vys <- c(vys, nb$mean_vy) }
      }
      svx[i] <- mean(vxs); svy[i] <- mean(vys)
    }
    fl[, c("svx", "svy") := .(svx, svy)]

    # Vorticity
    vort <- rep(NA_real_, nrow(fl))
    for (i in seq_len(nrow(fl))) {
      xb <- fl$x_bin[i]; yb <- fl$y_bin[i]
      e <- fl[.(xb + bs, yb)]; w <- fl[.(xb - bs, yb)]
      n_nb <- fl[.(xb, yb - bs)]; s_nb <- fl[.(xb, yb + bs)]
      dvy_dx <- if (nrow(e) == 1 && nrow(w) == 1)
                  (e$svy - w$svy) / (2 * bs) else NA_real_
      dvx_dy <- if (nrow(n_nb) == 1 && nrow(s_nb) == 1)
                  (s_nb$svx - n_nb$svx) / (2 * bs) else NA_real_
      if (!is.na(dvy_dx) && !is.na(dvx_dy)) vort[i] <- dvy_dx - dvx_dy
    }
    fl[, vorticity := vort]
    fl[, speed_2d := sqrt(svx^2 + svy^2)]
    fl[, swirl := abs(vorticity) * bs / pmax(speed_2d, 1e-8)]
    fl[, move_type := fifelse(
      !is.na(swirl) & swirl > swirl_thr, "Circular",
      fifelse(mean_vy > 0, "Downward (VP)", "Upward (AP)"))]
    fl[, move_type := factor(move_type,
        levels = c("Upward (AP)", "Circular", "Downward (VP)"))]
    fl[, species := sp_label]
    fl[, hour := h]
    res[[length(res) + 1]] <- fl
  }
  rbindlist(res, fill = TRUE)
}

cat("  Building hourly flow fields (this is the slowest step)...\n")
flow_m <- build_hourly_flow(vel_m, "Medaka")
flow_z <- build_hourly_flow(vel_z, "Zebrafish")
flow_all <- rbind(flow_m, flow_z, fill = TRUE)

# Hour facet labels — include n_bins and majority type
hour_majority <- flow_all[, {
  n <- .N
  maj <- names(sort(table(move_type), decreasing = TRUE))[1]
  pct <- max(table(move_type)) / n * 100
  .(majority = maj, maj_pct = pct, n_bins = n)
}, by = .(species, hour)]

cat("  Hourly majority movement type:\n")
print(hour_majority)
fwrite(hour_majority, file.path(OUT_DIR, "hourly_majority_movement.csv"))
fwrite(flow_all, file.path(OUT_DIR, "hourly_flow_fields.csv"))

# Pre-compute facet labels
flow_all <- merge(flow_all, hour_majority, by = c("species", "hour"))
flow_all[, hour_lab := sprintf("h%02d  maj: %s (%.0f%%)",
                                hour, majority, maj_pct)]

plot_flow_panel <- function(dt_sp) {
  sp <- unique(dt_sp$species)
  ggplot(dt_sp) +
    geom_segment(aes(x = x_bin, y = y_bin,
                     xend = x_bin + mean_vx * ARROW_SCALE,
                     yend = y_bin + mean_vy * ARROW_SCALE,
                     color = move_type),
                 arrow = arrow(length = unit(ARROW_HEAD, "cm"),
                               type = "closed"),
                 linewidth = ARROW_LW, alpha = 0.85) +
    facet_wrap(~ hour_lab, ncol = 4) +
    scale_color_manual(values = flow_cols, drop = FALSE, name = NULL) +
    scale_y_reverse() + coord_fixed() +
    labs(title = sprintf("%s — hourly XY flow fields", sp),
         subtitle = sprintf("bin = %d µm, min %d pts/bin | colour = bin movement type",
                            FLOW_BIN_UM, FLOW_MIN_N),
         x = "X (µm)", y = "Y (µm)  [AP top, VP bottom]") +
    theme_pub(10) +
    theme(strip.text = element_text(size = 8))
}

fig4_m <- plot_flow_panel(flow_all[species == "Medaka"])
fig4_z <- plot_flow_panel(flow_all[species == "Zebrafish"])

n_h_m <- uniqueN(flow_all[species == "Medaka", hour])
n_h_z <- uniqueN(flow_all[species == "Zebrafish", hour])
save_pdf(fig4_m, "04a_flow_per_hour_medaka.pdf",
         w = 18, h = max(6, 4 * ceiling(n_h_m / 4)))
save_pdf(fig4_z, "04b_flow_per_hour_zebrafish.pdf",
         w = 18, h = max(6, 4 * ceiling(n_h_z / 4)))

# Combined majority-type strip
maj_strip <- ggplot(hour_majority,
                    aes(hour, factor(species, levels = c("Zebrafish", "Medaka")),
                        fill = majority)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.0f%%", maj_pct)), size = 3) +
  scale_fill_manual(values = flow_cols, drop = FALSE, name = NULL) +
  labs(title = "Majority bin movement type per hour",
       x = "Hour (from start of acquisition)", y = NULL) +
  theme_pub() + theme(panel.grid = element_blank())
save_pdf(maj_strip, "04c_majority_movement_per_hour.pdf", w = 14, h = 4)

# =============================================================================
# FIGURE 5: MATCHED-TRACK COMPARISON (FAIR SAMPLING)
# =============================================================================
# Filter: 10 <= duration <= 60 min  AND  net displacement >= 20 µm.
# To compare two species sampled at different rates (medaka 30 s, zebrafish
# 120 s), we ALWAYS use a common 120-s, non-overlapping step interval here.
# Net displacement is from the endpoints of the full track (sampling-
# invariant); path length, turning, and MSD all use the 120-s steps so
# neither species accumulates extra jitter.  This gives one definitive
# answer for straightness, mean turning angle, and MSD.
# =============================================================================

banner("FIGURE 5 — matched-track comparison (fair 120-s sampling)")

# Endpoint-based per-track summary (sampling-invariant)
endpoint_per_track <- function(sp, fi_sec, sp_label) {
  fi_min <- fi_sec / 60
  setkey(sp, TRACK_ID, FRAME)
  sp[, .N, by = TRACK_ID][N >= 3, TRACK_ID] -> keep
  sp[TRACK_ID %in% keep,
     .(n_full       = .N,
       duration_min = (max(FRAME) - min(FRAME)) * fi_min,
       net_disp_um  = sqrt((POSITION_X[.N] - POSITION_X[1])^2 +
                           (POSITION_Y[.N] - POSITION_Y[1])^2 +
                           (POSITION_Z[.N] - POSITION_Z[1])^2)),
     by = TRACK_ID][, species := sp_label][]
}
ep_m <- endpoint_per_track(sp_m, MEDAKA_FI, "Medaka")
ep_z <- endpoint_per_track(sp_z, ZEB_FI,    "Zebrafish")
ep_all <- rbind(ep_m, ep_z)

matched_endpoint <- ep_all[duration_min >= 10 & duration_min <= 60 &
                              net_disp_um  >= 20]
ids_m <- matched_endpoint[species == "Medaka",    TRACK_ID]
ids_z <- matched_endpoint[species == "Zebrafish", TRACK_ID]
cat(sprintf("  Endpoint-matched tracks: Medaka n=%d, Zebrafish n=%d\n",
            length(ids_m), length(ids_z)))

# Fair-sampling: keep only frames at multiples of `retain_mod` so both species
# share a 120-s sampling grid.  Then derive BOTH path length and endpoint
# net displacement from the same sub-sampled positions, so straightness is
# always <= 1.
fair_subsample <- function(df, retain_mod) {
  d <- copy(df)[FRAME %% retain_mod == 0]
  setkey(d, TRACK_ID, FRAME)
  d
}
fair_per_track <- function(sub, fi_sec, lag, sp_label) {
  fi_min <- fi_sec / 60
  # step vectors between consecutive retained frames
  sub[, `:=`(
    dx = POSITION_X - shift(POSITION_X, 1),
    dy = POSITION_Y - shift(POSITION_Y, 1),
    dz = POSITION_Z - shift(POSITION_Z, 1),
    df_ = FRAME      - shift(FRAME,      1)
  ), by = TRACK_ID]
  sub[, ok := !is.na(df_) & df_ == lag]
  sub[, disp_3d := ifelse(ok, sqrt(dx^2 + dy^2 + dz^2), NA_real_)]
  # turning angle between consecutive valid steps
  sub[, `:=`(dx_p = shift(dx, 1), dy_p = shift(dy, 1), dz_p = shift(dz, 1),
              ok_p = shift(ok, 1)), by = TRACK_ID]
  sub[, mag_p := sqrt(dx_p^2 + dy_p^2 + dz_p^2)]
  sub[, dot_p := dx*dx_p + dy*dy_p + dz*dz_p]
  sub[ok & ok_p, turning_angle :=
        acos(pmin(pmax(dot_p / (disp_3d * mag_p), -1), 1)) * 180 / pi]
  # per-track summary on the SAME retained frames
  out <- sub[, {
    valid <- which(ok)
    if (length(valid) < 3) NULL else {
      # first and last retained positions (the ends of the contiguous
      # sub-sampled chain)
      idx_first <- min(c(valid, valid - 1L))   # frame before first valid step
      idx_last  <- max(valid)
      ndisp <- sqrt((POSITION_X[idx_last] - POSITION_X[idx_first])^2 +
                    (POSITION_Y[idx_last] - POSITION_Y[idx_first])^2 +
                    (POSITION_Z[idx_last] - POSITION_Z[idx_first])^2)
      tot   <- sum(disp_3d[valid])
      .(n_fair_steps   = length(valid),
        duration_min   = (FRAME[idx_last] - FRAME[idx_first]) * (fi_sec / 60),
        net_disp_um    = ndisp,
        total_path_um  = tot,
        straightness   = ndisp / tot,
        mean_turning_deg = mean(turning_angle, na.rm = TRUE))
    }
  }, by = TRACK_ID]
  out[, species := sp_label]
  out
}
sub_m <- fair_subsample(sp_m[TRACK_ID %in% ids_m], 4L)
sub_z <- fair_subsample(sp_z[TRACK_ID %in% ids_z], 1L)
matched <- rbind(
  fair_per_track(sub_m, MEDAKA_FI, lag = 4L, "Medaka"),
  fair_per_track(sub_z, ZEB_FI,    lag = 1L, "Zebrafish"))
matched <- matched[!is.na(straightness) & straightness <= 1.0001]
matched[, straightness := pmin(straightness, 1)]
matched[, species := factor(species, levels = c("Medaka", "Zebrafish"))]
fwrite(matched, file.path(OUT_DIR, "matched_tracks.csv"))

cat("  Per-track metrics (fair 120-s sampling):\n")
print(matched[, .(n              = .N,
                   med_dur        = median(duration_min),
                   med_disp       = median(net_disp_um),
                   med_straight   = median(straightness),
                   mean_straight  = mean(straightness),
                   med_turn_deg   = median(mean_turning_deg, na.rm = TRUE),
                   mean_turn_deg  = mean(mean_turning_deg, na.rm = TRUE)),
              by = species])

ws_s <- t.test(straightness     ~ species, data = matched)
ws_t <- t.test(mean_turning_deg ~ species, data = matched)
cat(sprintf("  Welch t  straightness: t=%+.2f, p=%.2e\n",
            ws_s$statistic, ws_s$p.value))
cat(sprintf("  Welch t  mean turning: t=%+.2f, p=%.2e\n",
            ws_t$statistic, ws_t$p.value))

p5a <- ggplot(matched, aes(species, straightness, fill = species)) +
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.18, outlier.size = 0.4, alpha = 0.8) +
  scale_fill_manual(values = species_colors) +
  labs(title = "Straightness",
       subtitle = sprintf("Welch p = %.2g  (fair 120-s steps)", ws_s$p.value),
       x = NULL, y = "net / total path") +
  theme_pub() + theme(legend.position = "none")

p5b <- ggplot(matched[!is.na(mean_turning_deg)],
              aes(species, mean_turning_deg, fill = species)) +
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.18, outlier.size = 0.4, alpha = 0.8) +
  scale_fill_manual(values = species_colors) +
  labs(title = "Mean turning angle",
       subtitle = sprintf("Welch p = %.2g  (fair 120-s steps)", ws_t$p.value),
       x = NULL, y = "deg between consecutive 120-s steps") +
  theme_pub() + theme(legend.position = "none")

# MSD at COMMON real-time lags (2, 4, 6, ... 60 min) for BOTH species.
# For medaka we accept only frame pairs whose Δframe is a multiple of 4
# (≡ 120 s), so both species are compared at identical real-time lags.
msd_matched_fair <- function(dt_sp, ids, fi_sec, frames_per_120s, max_k) {
  d <- dt_sp[TRACK_ID %in% ids]
  setkey(d, TRACK_ID, FRAME)
  out <- vector("list", max_k)
  for (k in seq_len(max_k)) {
    lag_f <- k * frames_per_120s
    pairs <- d[, {
      n <- .N
      if (n > lag_f) {
        i1 <- 1:(n - lag_f); i2 <- (1 + lag_f):n
        keep <- (FRAME[i2] - FRAME[i1]) == lag_f
        if (sum(keep) > 0) {
          dsq <- (POSITION_X[i2[keep]] - POSITION_X[i1[keep]])^2 +
                 (POSITION_Y[i2[keep]] - POSITION_Y[i1[keep]])^2 +
                 (POSITION_Z[i2[keep]] - POSITION_Z[i1[keep]])^2
          .(msd = mean(dsq), n = length(dsq))
        } else .(msd = numeric(0), n = integer(0))
      } else .(msd = numeric(0), n = integer(0))
    }, by = TRACK_ID]
    if (nrow(pairs) > 0 && sum(pairs$n) > 0) {
      out[[k]] <- data.table(lag_min  = k * 2,
                              mean_msd = weighted.mean(pairs$msd, pairs$n),
                              n_pairs  = sum(pairs$n))
    }
  }
  rbindlist(out)
}
max_k <- 30L   # 60 min
msd_m <- msd_matched_fair(sp_m, ids_m, MEDAKA_FI, 4L, max_k)[, species := "Medaka"]
msd_z <- msd_matched_fair(sp_z, ids_z, ZEB_FI,    1L, max_k)[, species := "Zebrafish"]
msd <- rbind(msd_m, msd_z)
fwrite(msd, file.path(OUT_DIR, "matched_msd.csv"))

# Power-law slope (alpha) of MSD vs lag, 2..30 min: alpha=2 ballistic, alpha=1 diffusive
msd_fit <- msd[lag_min >= 2 & lag_min <= 30,
                .(alpha = coef(lm(log(mean_msd) ~ log(lag_min)))[2]),
                by = species]
cat("  MSD power-law alpha (lag 2-30 min):\n"); print(msd_fit)

p5c <- ggplot(msd[lag_min <= 60],
              aes(lag_min, mean_msd, color = species)) +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  scale_color_manual(values = species_colors) +
  scale_x_log10() + scale_y_log10(labels = label_comma()) +
  geom_text(data = msd_fit,
            aes(x = 3, y = c(8e2, 3e3)[match(species, c("Medaka","Zebrafish"))],
                label = sprintf("alpha = %.2f", alpha),
                color = species),
            inherit.aes = FALSE, hjust = 0, size = 4, show.legend = FALSE) +
  labs(title = "MSD vs lag (log-log)",
       subtitle = sprintf("matched tracks (Medaka n=%d, Zebrafish n=%d) | common 120-s lags",
                          length(ids_m), length(ids_z)),
       x = "lag (min)", y = expression("MSD (" * mu * "m" ^ 2 * ")"),
       color = NULL) +
  theme_pub()

p5d <- ggplot(msd[lag_min <= 60],
              aes(lag_min, mean_msd, color = species)) +
  geom_line(linewidth = 1) + geom_point(size = 1.6) +
  scale_color_manual(values = species_colors) +
  labs(title = "MSD vs lag (linear)",
       x = "lag (min)", y = expression("MSD (" * mu * "m" ^ 2 * ")"),
       color = NULL) +
  theme_pub()

fig5 <- (p5a + p5b) / (p5c + p5d) +
  plot_annotation(
    title = "Matched-track comparison (fair 120-s sampling)",
    subtitle = "Tracks: 10-60 min duration AND >=20 um net displacement; path length, turning, MSD all at common 120-s steps",
    theme = theme(plot.title = element_text(face = "bold", size = 14)))
save_pdf(fig5, "05_matched_tracks.pdf", w = 13, h = 11)

banner("DONE")
cat(sprintf("  outputs in %s/\n", OUT_DIR))
