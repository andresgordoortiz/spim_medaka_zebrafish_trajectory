# =============================================================================
# ONE-SHOT MEDAKA vs ZEBRAFISH COMPARISON
# =============================================================================
# Generates the plot set requested:
#   (1) head-mesoderm / bulge cell density over time (medaka/zebrafish,
#       aligned to ingression)
#   (2) movement direction over time (epiboly fraction, mean v_epiboly,
#       per-step type composition)
#   (3) thickness in same plot (both species)
#   (4) flow fields per hour, both species, + majority movement type per hour
#   (5) track comparison (20-60 min duration, >=20 um net displacement):
#       straightness, mean turning angle, MSD + HTML QC viewer
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

MEDAKA_INGRESSION_FRAME <- 199L
ZEB_INGRESSION_FRAME    <- 40L
FIG1_ZEB_TCAP_MIN       <- 170

MEDAKA_VOXEL_UM <- 1.05152
ZEB_VOXEL_UM    <- 1.24785

MEDAKA_MIN_FRAMES <- 20L   # 10 min
ZEB_MIN_FRAMES    <- 5L    # 10 min

MEDAKA_VEL_LAG    <- 4L    # 2 min
ZEB_VEL_LAG       <- 1L
MEDAKA_SMOOTH_K   <- 5L
ZEB_SMOOTH_K      <- 3L

MARGIN_WIDTH_DEG  <- 5     # band ±5° around margin landmark

FIG1_ZEB_STRIP_MIN_START_CELLS   <- 200L
FIG1_ZEB_STRIP_MIN_STABLE_RATIO  <- 0.90
FIG1_ZEB_STRIP_SEARCH_PAD_DEG    <- 8L
FIG1_COMPARE_BIN_MIN             <- 2

FIG5_MIN_DURATION_MIN            <- 20
FIG5_MAX_DURATION_MIN            <- 60
FIG5_MIN_NET_DISP_UM             <- 20
FIG5_QC_SAMPLE_TRACKS            <- 250L
FIG5_QC_SEED                     <- 1L
FIG5_QC_TRAIL_STEPS              <- 6L

FLOW_BIN_UM       <- 30
FLOW_MIN_N        <- 10
ARROW_SCALE       <- 60
ARROW_HEAD        <- 0.14
ARROW_LW          <- 0.7
VORT_SWIRL_THRESH <- 0.3

species_colors <- c("Medaka" = "#E69F00", "Zebrafish" = "#0072B2")
step_type_cols <- c("Epiboly"     = "#D95A4E",
                    "Convergence" = "#4B9A68",
                    "Animalward"  = "#4C7FB8",
                    "Ingression"  = "#8A5FA8")
STEP_TYPE_LEVELS <- c("Epiboly", "Convergence", "Animalward", "Ingression")
flow_cols <- c("Upward (AP)" = "#2166AC", "Circular" = "#1A9850",
               "Downward (VP)" = "#B2182B")

FIG1_TIME_META <- data.table(
  species = c("Medaka", "Zebrafish"),
  ingression_frame = c(MEDAKA_INGRESSION_FRAME, ZEB_INGRESSION_FRAME),
  fi_sec = c(MEDAKA_FI, ZEB_FI)
)
FIG1_TIME_META[, ingression_time_min := ingression_frame * fi_sec / 60]
FIG1_BASELINE_FROM_INGR_MIN <- -30

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
  tmp <- tempfile(pattern = "plot_", fileext = ".pdf")
  on.exit(unlink(tmp), add = TRUE)
  ggsave(tmp, p, width = w, height = h, device = "pdf")
  if (!file.exists(tmp)) {
    stop("pdf device did not write output: ", tmp)
  }
  ok <- file.copy(tmp, path, overwrite = TRUE)
  if (!ok) {
    stop("failed to overwrite output pdf: ", path)
  }
  cat(sprintf("  saved %s\n", path))
}

save_html_file <- function(html_lines, name) {
  path <- file.path(OUT_DIR, name)
  writeLines(enc2utf8(html_lines), path, useBytes = TRUE)
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

read_landmark <- function(dir, name) {
  d <- fread(file.path(dir, "gastrulation_landmarks.csv"))
  row <- d[landmark == name]
  if (nrow(row) == 0) return(c(theta = NA_real_, phi = NA_real_))
  c(theta = as.numeric(row$theta[1]), phi = as.numeric(row$phi[1]))
}
lm_m_margin <- read_landmark(MEDAKA_INPUT, "margin")
lm_z_margin <- read_landmark(ZEBRAFISH_INPUT, "margin")
MARGIN_M <- lm_m_margin["theta"]
MARGIN_Z <- lm_z_margin["theta"]
lm_m_ingr <- read_landmark(MEDAKA_INPUT, "ingression_center")
lm_z_ingr <- read_landmark(ZEBRAFISH_INPUT, "ingression_center")
cat(sprintf("  Margin θ:           medaka %.2f°, zebrafish %.2f°\n", MARGIN_M, MARGIN_Z))

# Data-driven bulge / ingression-cluster detection (same algorithm as
# gastrulation_dynamics_medaka.R step 2): tracks whose max depth is in the top
# 5 % live in the cluster; centre = mean (θ,φ) of those tracks; radius =
# 2·sd(angular distance) + 2°, capped at 15°.  Zebrafish "shield" and medaka
# "ingression bulge" are biologically the same structure (site of cell
# internalisation at the dorsal margin).
INGR_DEPTH_PCTILE <- 0.95
detect_bulge <- function(sp, fallback_theta, fallback_phi) {
  # Where does the deep tissue ACCUMULATE late in the recording?  Use only
  # frames from the second half of the recording (after t_max/2) and pick the
  # deepest 5% of cells in that window; the centre is the spatial mean of
  # those observations.  This avoids the bias from tracks that dive early and
  # then drift away \u2014 we want the position where the bulge actually sits at
  # the end of the recording.
  t_max <- max(sp$FRAME, na.rm = TRUE)
  late  <- sp[FRAME >= t_max / 2 & is.finite(SPHERICAL_DEPTH)]
  thr <- as.numeric(quantile(late$SPHERICAL_DEPTH, INGR_DEPTH_PCTILE,
                              na.rm = TRUE))
  cand <- late[SPHERICAL_DEPTH >= thr]
  if (nrow(cand) >= 50) {
    ct <- mean(cand$THETA_DEG, na.rm = TRUE)
    cp <- mean(cand$PHI_DEG,   na.rm = TRUE)
    cand[, ang_dist := sqrt((THETA_DEG - ct)^2 + (PHI_DEG - cp)^2)]
    # Tight disc: 75th percentile of angular distance from the centre, capped
    # at 7 deg.  Earlier we used 2*sd + 2 deg (cap 15) which produced a disc
    # large enough to engulf normal-density tissue around the deep core and
    # diluted the bulge-density signal.
    rad <- min(as.numeric(quantile(cand$ang_dist, 0.75, na.rm = TRUE)), 7)
    list(theta = ct, phi = cp, radius = rad, n_cand = nrow(cand),
         depth_thresh = thr,
         source = "data-driven (deepest 5% in late half of recording)")
  } else {
    list(theta = fallback_theta, phi = fallback_phi, radius = 10,
         n_cand = nrow(cand), depth_thresh = thr,
         source = "landmark fallback")
  }
}
BULGE_M <- detect_bulge(sp_m, lm_m_ingr["theta"], lm_m_ingr["phi"])
BULGE_Z <- detect_bulge(sp_z, lm_z_ingr["theta"], lm_z_ingr["phi"])
cat(sprintf("  Bulge medaka:    θ=%.1f° φ=%.1f° radius=%.1f°  (%d deep tracks, depth ≥ %.1f µm, %s)\n",
            BULGE_M$theta, BULGE_M$phi, BULGE_M$radius,
            BULGE_M$n_cand, BULGE_M$depth_thresh, BULGE_M$source))
cat(sprintf("  Bulge zebrafish: θ=%.1f° φ=%.1f° radius=%.1f°  (%d deep tracks, depth ≥ %.1f µm, %s)\n",
            BULGE_Z$theta, BULGE_Z$phi, BULGE_Z$radius,
            BULGE_Z$n_cand, BULGE_Z$depth_thresh, BULGE_Z$source))

select_stable_zeb_strip <- function(sp, bulge_theta, landmark_theta,
                                    width = MARGIN_WIDTH_DEG,
                                    ing_frame = ZEB_INGRESSION_FRAME,
                                    min_start_cells = FIG1_ZEB_STRIP_MIN_START_CELLS,
                                    min_stable_ratio = FIG1_ZEB_STRIP_MIN_STABLE_RATIO,
                                    search_pad_deg = FIG1_ZEB_STRIP_SEARCH_PAD_DEG) {
  start_frame <- min(sp$FRAME, na.rm = TRUE)
  centers <- seq(floor(bulge_theta) - search_pad_deg,
                 floor(landmark_theta), by = 1)
  cand <- rbindlist(lapply(centers, function(center) {
    lo <- center - width
    hi <- center + width
    by_frame <- sp[, .(
      n_total = .N,
      n_band = sum(THETA_DEG >= lo & THETA_DEG <= hi)
    ), by = FRAME]
    start_row <- by_frame[FRAME == start_frame]
    ing_row   <- by_frame[FRAME == ing_frame]
    if (!nrow(start_row) || !nrow(ing_row) || ing_row$n_band[1] == 0) {
      return(NULL)
    }
    frac_start <- start_row$n_band[1] / start_row$n_total[1]
    frac_ing   <- ing_row$n_band[1]   / ing_row$n_total[1]
    data.table(
      center = center,
      n_start = start_row$n_band[1],
      n_ing = ing_row$n_band[1],
      frac_start = frac_start,
      frac_ing = frac_ing,
      stable_ratio = frac_start / frac_ing
    )
  }), fill = TRUE)

  keep <- cand[n_start >= min_start_cells & stable_ratio >= min_stable_ratio]
  if (nrow(keep)) {
    best <- keep[which.max(center)]
    best[, selection_rule := sprintf(
      "most vegetal stable center (start frac >= %.0f%% of ingression frac; start n >= %d)",
      min_stable_ratio * 100, min_start_cells)]
  } else {
    best <- cand[which.max(frac_start)]
    best[, selection_rule := "fallback: highest starting occupancy fraction"]
  }
  best[, width := width]
  best[]
}

FIG1_ZEB_STRIP <- select_stable_zeb_strip(sp_z, BULGE_Z$theta, MARGIN_Z)
FIG1_ZEB_STRIP_CENTER <- FIG1_ZEB_STRIP$center[1]
cat(sprintf(
  "  Figure 1 zebrafish comparison strip: θ=%.1f° ± %d°  (start %d cells / %.1f%% of field, ingression %d cells / %.1f%%, ratio %.3f; %s)\n",
  FIG1_ZEB_STRIP_CENTER, MARGIN_WIDTH_DEG,
  FIG1_ZEB_STRIP$n_start[1], 100 * FIG1_ZEB_STRIP$frac_start[1],
  FIG1_ZEB_STRIP$n_ing[1], 100 * FIG1_ZEB_STRIP$frac_ing[1],
  FIG1_ZEB_STRIP$stable_ratio[1], FIG1_ZEB_STRIP$selection_rule[1]))

# Bulge disc outlines (used by both Fig 1 and Fig 3 heatmap overlays)
bulge_circle_poly <- function(b, species_label, n = 80) {
  ang <- seq(0, 2 * pi, length.out = n)
  data.table(species = factor(species_label, levels = c("Medaka", "Zebrafish")),
             phi   = b$phi   + b$radius * cos(ang),
             theta = b$theta + b$radius * sin(ang))
}
bulge_poly <- rbind(bulge_circle_poly(BULGE_M, "Medaka"),
                    bulge_circle_poly(BULGE_Z, "Zebrafish"))

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
                "RADIAL_DIST", "THETA_DEG", "PHI_DEG",
                "SPHERICAL_DEPTH")) {
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
  dt[, dDepth := SPHERICAL_DEPTH_SM - shift(SPHERICAL_DEPTH_SM, vel_lag),
     by = TRACK_ID]
  dt[, `:=`(
    inst_speed    = disp_3d / (dt_f * fi_min),
    vx_um_min     = dx / (dt_f * fi_min),
    vy_um_min     = dy / (dt_f * fi_min),
    v_epiboly     = RADIAL_DIST * (dTheta * pi / 180) / (dt_f * fi_min),
    v_convergence = RADIAL_DIST * sin(theta_rad) *
                    (dPhi * pi / 180) / (dt_f * fi_min),
    v_depth       = dDepth / (dt_f * fi_min)  # positive = deepening (inward)
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

# Per-step direction classification (3 categories — Convergence / Epiboly /
# Animalward).  Ingression cannot be defined per step robustly (single inward
# steps are noisy); we use the per-track classifier below for the composition
# stack so it matches the previous gastrulation_dynamics_05 plot.
classify_steps <- function(dt) {
  dt[, step_type := fcase(
    abs(v_convergence) > abs(v_epiboly), "Convergence",
    v_epiboly > 0,                       "Epiboly",
    default                            = "Animalward"
  )]
  dt[, step_type := factor(step_type,
                            levels = c("Epiboly", "Convergence", "Animalward"))]
  dt[]
}

vel_m <- classify_steps(vel_m)
vel_z <- classify_steps(vel_z)
vel_m[, species := "Medaka"];  vel_z[, species := "Zebrafish"]
vel_m[, time_min := FRAME * MEDAKA_FI / 60]
vel_z[, time_min := FRAME * ZEB_FI    / 60]

cat(sprintf("  velocity steps: medaka %s, zebrafish %s\n",
            format(nrow(vel_m), big.mark = ","),
            format(nrow(vel_z), big.mark = ",")))

# -----------------------------------------------------------------------------
# Per-TRACK movement_type — verbatim port of gastrulation_dynamics_medaka.R
# (the classifier that produced medaka_dynamics_05_spatial_flow).  Two key
# pieces of physics:
#
#   * v_epiboly      = R · (Δθ/Δt) in µm/min along the meridian — positive
#                      means the cell moves vegetal-ward (toward the equator,
#                      i.e. the epiboly direction).
#   * v_convergence  = R · sin(θ) · (Δφ/Δt) in µm/min along the equator — its
#                      sign on its own is meaningless because cells on the two
#                      flanks converge from opposite sides.
#   * v_conv_dorsal  = −v_convergence · sign(φ − φ_bulge): positive whenever
#                      the step shortens the equatorial gap to the dorsal
#                      bulge / shield axis.  THIS is the convergence signal.
#   * v_depth        = Δ(spherical depth)/Δt; positive = moving inward.
#
# Categories (in priority order):
#   Ingression  : track ENDS in the ingression region AND has gained more depth
#                 than the global P90 of (last-first) depth change.
#                 Region differs by species — see ingr_region_* below.
#   Convergence : |mean_v_conv_dorsal| > |mean_v_epiboly|  AND  > 0.
#   Epiboly     : mean_v_epiboly > 0.
#   Animalward  : otherwise (rare; mostly noise / animal-pole drift).
# -----------------------------------------------------------------------------
INGR_DEPTH_CHANGE_PCTILE <- 0.90

# Ingression spatial region per species:
#  - Medaka  : a localized DISC = the data-driven bulge cluster.
#              Cells internalise at one spot (the germ-ring bulge).
#  - Zebrafish: the whole MARGIN STRIP (no localized bulge); cells involute
#              all around the ring.  The shield is the dorsal hot spot but
#              ingression is not restricted to it.
INGR_REGION_M <- list(type = "disc",  theta = BULGE_M$theta, phi = BULGE_M$phi,
                       radius = BULGE_M$radius)
# Half-width for the zebrafish margin strip used in ingression classification
# (slightly wider than the thickness-zone margin so we catch involution that
# already passed the strict band).
INGR_STRIP_HALFW_Z <- 12
INGR_REGION_Z <- list(type = "strip", margin_theta = MARGIN_Z,
                       half_width = INGR_STRIP_HALFW_Z)

in_ingr_region <- function(theta, phi, r) {
  if (r$type == "disc") {
    sqrt((theta - r$theta)^2 + (phi - r$phi)^2) < r$radius
  } else {
    abs(theta - r$margin_theta) <= r$half_width
  }
}

classify_tracks <- function(vel, region, dorsal_phi) {
  # per-step dorsal-directed convergence (positive = toward the dorsal axis)
  vel[, conv_dorsal_step := -v_convergence * sign(PHI_DEG - dorsal_phi)]
  tm <- vel[!is.na(v_epiboly), .(
    mean_v_epiboly     = mean(v_epiboly,        na.rm = TRUE),
    mean_conv_dorsal   = mean(conv_dorsal_step, na.rm = TRUE),
    mean_v_depth       = mean(v_depth,          na.rm = TRUE),
    end_theta          = last(THETA_DEG),
    end_phi            = last(PHI_DEG),
    depth_change       = last(SPHERICAL_DEPTH) - first(SPHERICAL_DEPTH),
    max_depth          = max(SPHERICAL_DEPTH,   na.rm = TRUE),
    n_steps            = .N
  ), by = TRACK_ID]
  thr_depth <- as.numeric(quantile(tm$depth_change,
                                    INGR_DEPTH_CHANGE_PCTILE, na.rm = TRUE))
  if (!is.finite(thr_depth) || thr_depth <= 0) thr_depth <- 5
  tm[, in_ingr_end := in_ingr_region(end_theta, end_phi, region)]
  tm[, movement_type := fcase(
    in_ingr_end & depth_change > thr_depth,                             "Ingression",
    abs(mean_conv_dorsal) > abs(mean_v_epiboly) & mean_conv_dorsal > 0, "Convergence",
    mean_v_epiboly > 0,                                                 "Epiboly",
    default                                                             = "Animalward"
  )]
  vel[, conv_dorsal_step := NULL]
  list(track_class = tm[, .(TRACK_ID, movement_type, depth_change,
                             mean_v_depth, in_ingr_end)],
       depth_thresh = thr_depth)
}

trk_m <- classify_tracks(vel_m, INGR_REGION_M, BULGE_M$phi)
trk_z <- classify_tracks(vel_z, INGR_REGION_Z, BULGE_Z$phi)
cat(sprintf("  Track-level ingression depth threshold: medaka %.1f µm, zebrafish %.1f µm\n",
            trk_m$depth_thresh, trk_z$depth_thresh))
cat("  Track counts (Medaka): \n");    print(trk_m$track_class[, .N, by = movement_type])
cat("  Track counts (Zebrafish): \n"); print(trk_z$track_class[, .N, by = movement_type])

vel_m <- merge(vel_m, trk_m$track_class[, .(TRACK_ID, movement_type)],
                by = "TRACK_ID", all.x = TRUE)
vel_z <- merge(vel_z, trk_z$track_class[, .(TRACK_ID, movement_type)],
                by = "TRACK_ID", all.x = TRUE)
vel_m[, movement_type := factor(movement_type, levels = STEP_TYPE_LEVELS)]
vel_z[, movement_type := factor(movement_type, levels = STEP_TYPE_LEVELS)]

# =============================================================================
# FIGURE 1: BULGE NUCLEI DENSITY  ---  is the bulge a high-density region?
# =============================================================================
#
# We answer ONE question: does the bulge column contain more nuclei (per unit
# embryo-surface area) than other regions, and does this rise over time?
#
# What "bulge disc" means here:
#   For each species we define a small disc on the EMBRYO SURFACE around the
#   data-driven bulge centre.  The density measurement counts ALL nuclei whose
#   (theta, phi) falls inside that disc, AT ANY DEPTH (= the full radial
#   column under the disc, from outer surface to inner edge of the tracked
#   tissue).  We then divide that count by the disc's surface area on the
#   sphere.  So the units are "nuclei per 1000 um^2 of embryo surface" and
#   the metric is depth-INTEGRATED -- if the column has more cells stacked
#   inward, the number goes up.
#
# To make this concrete, panel A shows the disc on the embryo surface, and
# panel B shows a CROSS-SECTION (theta vs depth, longitudinal slice through
# the bulge centre) so the radial column under the disc is visible.  Panel C
# is the time series.
# =============================================================================

banner("FIGURE 1 - head mesoderm / bulge nuclei density")

# ---- Margin-strip per-frame density (kept for the CSV record only) ---------
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
md_z <- margin_density(sp_z, R_Z, FIG1_ZEB_STRIP_CENTER, PHI_RNG_Z, ZEB_FI)
md_m[, species := "Medaka"];  md_z[, species := "Zebrafish"]

# ---- Areas ---------------------------------------------------------------
imaged_area_um2 <- function(dt, R) {
  th_q <- as.numeric(quantile(dt$THETA_DEG, c(0.02, 0.98), na.rm = TRUE)) * pi / 180
  ph_q <- as.numeric(quantile(dt$PHI_DEG,   c(0.02, 0.98), na.rm = TRUE)) * pi / 180
  R^2 * (th_q[2] - th_q[1]) * sin(mean(th_q)) * (ph_q[2] - ph_q[1])
}
A_GLOBAL_M <- imaged_area_um2(sp_m, R_M)
A_GLOBAL_Z <- imaged_area_um2(sp_z, R_Z)
# Disc area: the in_disc selector uses Euclidean (theta,phi) distance, so the
# selected region is an ellipse on the sphere with semi-axes R*r_rad (theta)
# and R*sin(theta0)*r_rad (phi).  Area = pi * R^2 * r_rad^2 * sin(theta0).
# This is consistent with the selector (NOT the spherical-cap area, which
# would be pi*R^2*r_rad^2 -- 0.5-2% larger here).
ellipse_disc_area <- function(R, b)
  pi * R^2 * (b$radius * pi / 180)^2 * sin(b$theta * pi / 180)
A_BULGE_M  <- ellipse_disc_area(R_M, BULGE_M)
A_BULGE_Z  <- ellipse_disc_area(R_Z, BULGE_Z)
A_STRIP_M  <- md_m$area_um2[1]
A_STRIP_Z  <- md_z$area_um2[1]

# Control disc -- SAME LATITUDE (theta) as the bulge so we are NOT comparing
# the marginal band against the animal cap (which would be unfair: the imaged
# region barely reaches theta<65 deg).  We slide along phi to find an azimuth
# inside the imaged window, at least 2*radius away from the bulge, that
# carries the largest cell count over time at the bulge's theta -- that is,
# the most data-rich non-bulge location on the same latitude.  Same disc
# radius, so identical spherical-cap area, so densities are directly
# comparable nuclei-per-unit-area numbers.
pick_ctrl_phi <- function(sp, b, sep = NULL) {
  if (is.null(sep)) sep <- 2 * b$radius
  ring <- sp[!is.na(PHI_DEG) & abs(THETA_DEG - b$theta) < b$radius]
  if (!nrow(ring)) return(b$phi + 60)
  phi_lo <- as.numeric(quantile(ring$PHI_DEG, 0.02, na.rm = TRUE)) + b$radius
  phi_hi <- as.numeric(quantile(ring$PHI_DEG, 0.98, na.rm = TRUE)) - b$radius
  cand_phi <- seq(phi_lo, phi_hi, by = 2)
  cand_phi <- cand_phi[abs(cand_phi - b$phi) >= sep]
  if (!length(cand_phi)) return(b$phi + 60)
  counts <- vapply(cand_phi, function(p) sum(
    sqrt((ring$THETA_DEG - b$theta)^2 + (ring$PHI_DEG - p)^2) < b$radius),
    integer(1))
  cand_phi[which.max(counts)]
}
CTRL_M <- list(theta = BULGE_M$theta, phi = pick_ctrl_phi(sp_m, BULGE_M),
               radius = BULGE_M$radius)
CTRL_Z <- list(theta = BULGE_Z$theta, phi = pick_ctrl_phi(sp_z, BULGE_Z),
               radius = BULGE_Z$radius)
A_CTRL_M <- ellipse_disc_area(R_M, CTRL_M)   # = A_BULGE_M (same theta, same r)
A_CTRL_Z <- ellipse_disc_area(R_Z, CTRL_Z)   # = A_BULGE_Z (same theta, same r)

cat(sprintf("  Disc params -- Medaka:    bulge theta=%.1f phi=%.1f r=%.1f deg, A=%.1f x10^3 um^2 ; control theta=%.1f phi=%.1f\n",
            BULGE_M$theta, BULGE_M$phi, BULGE_M$radius, A_BULGE_M/1000, CTRL_M$theta, CTRL_M$phi))
cat(sprintf("  Disc params -- Zebrafish: bulge theta=%.1f phi=%.1f r=%.1f deg, A=%.1f x10^3 um^2 ; control theta=%.1f phi=%.1f\n",
            BULGE_Z$theta, BULGE_Z$phi, BULGE_Z$radius, A_BULGE_Z/1000, CTRL_Z$theta, CTRL_Z$phi))
cat(sprintf("  Margin strip area: medaka %.1f x10^3 um^2 ; zebrafish %.1f x10^3 um^2\n",
            A_STRIP_M/1000, A_STRIP_Z/1000))
cat(sprintf("  Imaged surface:    medaka %.1f x10^3 um^2 ; zebrafish %.1f x10^3 um^2\n",
            A_GLOBAL_M/1000, A_GLOBAL_Z/1000))

# ---- Selectors and helpers (time-cap independent) --------------------------
in_disc <- function(sp, b) sp[!is.na(PHI_DEG) &
  sqrt((THETA_DEG - b$theta)^2 + (PHI_DEG - b$phi)^2) < b$radius]
in_strip <- function(sp, m) sp[THETA_DEG >= m - MARGIN_WIDTH_DEG &
                                THETA_DEG <= m + MARGIN_WIDTH_DEG]
per_frame_density <- function(cells, A, fi, sp_lbl, reg_lbl) {
  cells[, .(n = .N), by = FRAME][,
        .(species = sp_lbl, region = reg_lbl, FRAME,
          time_min = FRAME * fi / 60,
          n_cells = n,
          area_um2 = A,
          density_per_1000um2 = n / A * 1000)]
}
REGION_LVL <- c("Bulge disc", "Comparison strip", "Global (whole imaged surface)")
REGION_COL <- c("Bulge disc"                    = "#E31A1C",
                "Comparison strip"              = "#5AB4AC",
                "Global (whole imaged surface)" = "grey40")
SNAP_NMAX <- 30000
PHI_SLICE_W <- 8
sample_n <- function(d, n) if (nrow(d) > n) d[sample(.N, n)] else d
late_snap <- function(sp, fi, t_cap, w = 10) {
  d <- sp[FRAME * fi / 60 <= t_cap]
  tmax <- max(d$FRAME) * fi / 60
  d[FRAME * fi / 60 >= tmax - w]
}
disc_outline <- function(b, sp_lbl, kind, n = 80) {
  ang <- seq(0, 2 * pi, length.out = n)
  data.table(species = sp_lbl, kind = kind,
             phi   = b$phi   + b$radius * cos(ang),
             theta = b$theta + b$radius * sin(ang))
}
peak_window_mean <- function(t, y, w = 30) {
  ord <- order(t); t <- t[ord]; y <- y[ord]
  y_smooth <- vapply(t, function(tc) mean(y[t >= tc - w/2 & t <= tc + w/2]),
                     numeric(1))
  i <- which.max(y_smooth)
  list(peak_mean = mean(y[t >= t[i] - w/2 & t <= t[i] + w/2]),
       peak_time = t[i])
}

# ---- Fig 1 builder.  Uses species-specific caps and aligns time to ingression.
build_fig1 <- function(fname,
                       time_caps = c(Medaka = Inf, Zebrafish = FIG1_ZEB_TCAP_MIN),
                       title_suffix = "") {
  cap_m <- unname(time_caps[["Medaka"]])
  cap_z <- unname(time_caps[["Zebrafish"]])
  cat(sprintf("\n  -- build_fig1(medaka=%s, zebrafish=%s) -> %s\n",
              if (is.finite(cap_m)) sprintf("%.0f min", cap_m) else "full",
              if (is.finite(cap_z)) sprintf("%.0f min", cap_z) else "full",
              fname))
  # Per-frame densities (filtered by species-specific caps)
  sp_m_t <- sp_m[time_min <= cap_m]
  sp_z_t <- sp_z[time_min <= cap_z]
  d_ts <- rbind(
    per_frame_density(in_disc(sp_m_t, BULGE_M), A_BULGE_M, MEDAKA_FI, "Medaka",    "Bulge disc"),
    per_frame_density(in_disc(sp_z_t, BULGE_Z), A_BULGE_Z, ZEB_FI,    "Zebrafish", "Bulge disc"),
    per_frame_density(in_strip(sp_m_t, MARGIN_M), A_STRIP_M, MEDAKA_FI, "Medaka",    "Comparison strip"),
    per_frame_density(in_strip(sp_z_t, FIG1_ZEB_STRIP_CENTER), A_STRIP_Z, ZEB_FI,    "Zebrafish", "Comparison strip"),
    per_frame_density(sp_m_t, A_GLOBAL_M, MEDAKA_FI, "Medaka",    "Global (whole imaged surface)"),
    per_frame_density(sp_z_t, A_GLOBAL_Z, ZEB_FI,    "Zebrafish", "Global (whole imaged surface)"))
  d_ts <- merge(d_ts,
                FIG1_TIME_META[, .(species, ingression_frame,
                                   ingression_time_min)],
                by = "species", all.x = TRUE)
  d_ts[, time_from_ingression_min := time_min - ingression_time_min]
  d_ts[, species := factor(species, levels = c("Medaka", "Zebrafish"))]
  d_ts[, region := factor(region, levels = REGION_LVL)]

  # QC: drop frames where global density < 50% of species median
  global_qc <- d_ts[region == "Global (whole imaged surface)",
                    .(thr = 0.5 * median(density_per_1000um2)), by = species]
  bad <- d_ts[region == "Global (whole imaged surface)"][
    global_qc, on = "species"][density_per_1000um2 < thr, .(species, FRAME)]
  if (nrow(bad)) {
    cat(sprintf("     dropped %d broken frames (global < 50%% of median)\n", nrow(bad)))
    d_ts <- d_ts[!bad, on = c("species", "FRAME")]
  }

  density_qc <- d_ts[, .(
    area_um2 = first(area_um2),
    max_formula_error = max(abs(density_per_1000um2 - n_cells / area_um2 * 1000),
                            na.rm = TRUE),
    median_n_cells = as.numeric(median(n_cells))
  ), by = .(species, region)]
  cat("     density QC (density = n_cells / area_um2 * 1000):\n")
  print(density_qc)

  # Panel A: top view, last 10 min of the (capped) window
  snap_m <- late_snap(sp_m, MEDAKA_FI, cap_m)
  snap_z <- late_snap(sp_z, ZEB_FI,    cap_z)
  tA_m <- max(snap_m$time_min); tA_z <- max(snap_z$time_min)
  snap <- rbind(
    sample_n(snap_m, SNAP_NMAX)[, .(species = "Medaka",    PHI_DEG, THETA_DEG, SPHERICAL_DEPTH)],
    sample_n(snap_z, SNAP_NMAX)[, .(species = "Zebrafish", PHI_DEG, THETA_DEG, SPHERICAL_DEPTH)])
  snap[, species := factor(species, levels = c("Medaka", "Zebrafish"))]
  outlines <- rbind(
    disc_outline(BULGE_M, "Medaka",    "Bulge"),
    disc_outline(BULGE_Z, "Zebrafish", "Bulge"))
  outlines[, species := factor(species, levels = c("Medaka", "Zebrafish"))]
  band_rect_p1 <- data.table(
    species = factor(c("Medaka", "Zebrafish"), levels = c("Medaka", "Zebrafish")),
    ymin    = c(MARGIN_M - MARGIN_WIDTH_DEG,
                FIG1_ZEB_STRIP_CENTER - MARGIN_WIDTH_DEG),
    ymax    = c(MARGIN_M + MARGIN_WIDTH_DEG,
                FIG1_ZEB_STRIP_CENTER + MARGIN_WIDTH_DEG))

  pA <- ggplot(snap, aes(PHI_DEG, THETA_DEG)) +
    geom_rect(data = band_rect_p1, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
          fill = REGION_COL[["Comparison strip"]], alpha = 0.18) +
    geom_point(aes(color = SPHERICAL_DEPTH), size = 0.25, alpha = 0.55) +
    geom_polygon(data = outlines[kind == "Bulge"], inherit.aes = FALSE,
                 aes(phi, theta), fill = NA, color = "#E31A1C", linewidth = 0.9) +
    facet_wrap(~ species, scales = "free") +
    scale_y_reverse() +
    scale_color_viridis_c(option = "inferno", name = expression("depth ("*mu*"m)")) +
    labs(title = sprintf("A. Top view of the embryo surface (snapshot: medaka %.0f-%.0f min, zebrafish %.0f-%.0f min)",
                          tA_m - 10, tA_m, tA_z - 10, tA_z),
        subtitle = sprintf("Dot colour = depth (yellow = deep).  Red = head-mesoderm / bulge disc.  Teal band = comparison strip (medaka landmark margin; zebrafish stable strip at theta %.0f°).",
                 FIG1_ZEB_STRIP_CENTER),
         x = expression(varphi*" (deg)"), y = expression(theta*" (deg, animal pole up)")) +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 9, color = "grey25"))

  # Panel B: longitudinal cross-section through bulge centre
  slice_m <- late_snap(sp_m[abs(PHI_DEG - BULGE_M$phi) <= PHI_SLICE_W], MEDAKA_FI, cap_m)
  slice_z <- late_snap(sp_z[abs(PHI_DEG - BULGE_Z$phi) <= PHI_SLICE_W], ZEB_FI,    cap_z)
  classify_xs <- function(d, bulge) {
    d[, ang_b := sqrt((THETA_DEG - bulge$theta)^2 + (PHI_DEG - bulge$phi)^2)]
    d[, xs_region := fifelse(ang_b < bulge$radius, "Bulge column", "Elsewhere")]
  }
  classify_xs(slice_m, BULGE_M)
  classify_xs(slice_z, BULGE_Z)
  xs <- rbind(
    sample_n(slice_m, SNAP_NMAX)[, .(species = "Medaka",    THETA_DEG, SPHERICAL_DEPTH, xs_region)],
    sample_n(slice_z, SNAP_NMAX)[, .(species = "Zebrafish", THETA_DEG, SPHERICAL_DEPTH, xs_region)])
  xs[, species   := factor(species,   levels = c("Medaka", "Zebrafish"))]
  xs[, xs_region := factor(xs_region, levels = c("Elsewhere", "Bulge column"))]
  xs_bands <- rbind(
    data.table(species = "Medaka",    kind = "Bulge", xmin = BULGE_M$theta - BULGE_M$radius, xmax = BULGE_M$theta + BULGE_M$radius),
    data.table(species = "Zebrafish", kind = "Bulge", xmin = BULGE_Z$theta - BULGE_Z$radius, xmax = BULGE_Z$theta + BULGE_Z$radius))
  xs_bands[, species := factor(species, levels = c("Medaka", "Zebrafish"))]

  pB <- ggplot(xs, aes(THETA_DEG, SPHERICAL_DEPTH)) +
    geom_rect(data = xs_bands[kind == "Bulge"], inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "#E31A1C", alpha = 0.10) +
    geom_point(aes(color = xs_region), size = 0.25, alpha = 0.55) +
    scale_color_manual(values = c("Elsewhere" = "grey60",
                                  "Bulge column"   = "#E31A1C"), name = NULL) +
    facet_wrap(~ species, scales = "free") +
    scale_y_reverse() +
    labs(title = sprintf("B. Cross-section through bulge centre (slice |phi-phi_bulge| <= %g deg)", PHI_SLICE_W),
         subtitle = "Density counts every nucleus in the red column, all depths.  y-axis inverted: deeper = down.",
         x = expression(theta*" (deg)"), y = expression("spherical depth ("*mu*"m)")) +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 9, color = "grey25"))

    # Panel C: density vs real acquisition time
    pC <- ggplot(d_ts, aes(time_min, density_per_1000um2,
                         color = region)) +
    geom_point(alpha = 0.15, size = 0.4) +
    geom_smooth(method = "loess", span = 0.3, se = TRUE, linewidth = 1) +
    facet_wrap(~ species, nrow = 1, scales = "free_x") +
    scale_color_manual(values = REGION_COL, name = NULL) +
      labs(title = "C. Regional nuclei density trajectories over real time",
        subtitle = sprintf("Raw acquisition time in minutes.  Zebrafish comparison strip is the most vegetal band already in view before ingression (theta %.0f° ± %d°); panel D handles ingression alignment.",
                 FIG1_ZEB_STRIP_CENTER, MARGIN_WIDTH_DEG),
        x = "time (min)",
         y = expression("nuclei / 1000 "*mu*"m"^2)) +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 9, color = "grey25"),
          legend.position = "bottom")

  ratio_spec <- data.table(
    species = c("Medaka", "Zebrafish"),
      numerator_region = c("Bulge disc", "Bulge disc"),
      ratio_label = c("Head mesoderm disc / global",
              "Head mesoderm disc / global")
  )
  baseline_dt <- d_ts[, {
    dens <- density_per_1000um2[
      time_from_ingression_min >= FIG1_BASELINE_FROM_INGR_MIN &
      time_from_ingression_min < 0
    ]
    if (!length(dens)) dens <- density_per_1000um2[time_from_ingression_min < 0]
    if (!length(dens)) dens <- head(density_per_1000um2, 5)
    .(baseline_density = mean(dens, na.rm = TRUE))
  }, by = .(species, region)]
  numerator_dt <- merge(
    d_ts,
    ratio_spec,
    by.x = c("species", "region"),
    by.y = c("species", "numerator_region"),
    all = FALSE
  )
  numerator_dt <- merge(
    numerator_dt,
    baseline_dt,
    by = c("species", "region"),
    all.x = TRUE
  )
  setnames(numerator_dt,
           c("n_cells", "area_um2", "density_per_1000um2", "baseline_density"),
           c("numerator_n_cells", "numerator_area_um2",
             "numerator_density_per_1000um2", "numerator_baseline_density"))
  global_dt <- d_ts[region == "Global (whole imaged surface)",
                    .(species, FRAME, time_min, time_from_ingression_min,
                      global_n_cells = n_cells,
                      global_area_um2 = area_um2,
                      global_density_per_1000um2 = density_per_1000um2)]
  global_base <- baseline_dt[region == "Global (whole imaged surface)",
                             .(species, global_baseline_density = baseline_density)]
  ratio_ts <- merge(
    numerator_dt[, .(species, region, ratio_label, FRAME, time_min,
                     time_from_ingression_min,
                     numerator_n_cells, numerator_area_um2,
                     numerator_density_per_1000um2,
                     numerator_baseline_density)],
    global_dt,
    by = c("species", "FRAME", "time_min", "time_from_ingression_min"),
    all = FALSE
  )
  ratio_ts <- merge(ratio_ts, global_base, by = "species", all.x = TRUE)
  ratio_ts[, numerator_fold := numerator_density_per_1000um2 / numerator_baseline_density]
  ratio_ts[, global_fold := global_density_per_1000um2 / global_baseline_density]
  ratio_ts[, fold_change_over_global := numerator_fold / global_fold]

  ratio_plot_dt <- ratio_ts[, .(
    fold_change_over_global = mean(fold_change_over_global, na.rm = TRUE),
    numerator_density_per_1000um2 = mean(numerator_density_per_1000um2,
                                         na.rm = TRUE),
    global_density_per_1000um2 = mean(global_density_per_1000um2,
                                      na.rm = TRUE),
    n_obs = .N
  ), by = .(
    species,
    ratio_label,
    time_plot_min = floor(time_from_ingression_min / FIG1_COMPARE_BIN_MIN) *
      FIG1_COMPARE_BIN_MIN
  )]

  ratio_summary <- ratio_plot_dt[, .(
    baseline_window = sprintf("[%d, 0) min", FIG1_BASELINE_FROM_INGR_MIN),
    peak_ratio = round(max(fold_change_over_global, na.rm = TRUE), 3),
    peak_time_from_ingression = time_plot_min[which.max(fold_change_over_global)]
  ), by = .(species, ratio_label)]
  cat("     aligned fold-change / global summary:\n")
  print(ratio_summary)

  pD <- ggplot(ratio_plot_dt,
               aes(time_plot_min, fold_change_over_global,
                   color = species)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey65") +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey55") +
    geom_line(linewidth = 0.9, alpha = 0.95) +
    geom_point(size = 1.3, alpha = 0.85) +
    scale_color_manual(values = species_colors, name = NULL) +
    labs(title = "D. Same head-mesoderm disc in both species",
         subtitle = sprintf("Both curves use the same red disc definition per species, normalized as local fold change / global fold change.  2-min aligned bins remove the sampling-rate mismatch; values > 1 mean disc-specific enrichment beyond embryo-wide densification.  Baseline window = [%d, 0) min.",
                            FIG1_BASELINE_FROM_INGR_MIN),
         x = "time from ingression (min)",
         y = "head-mesoderm fold change / global fold change") +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 9, color = "grey25"))

  fig <- (pA / pB / pC / pD) +
    plot_layout(heights = c(1.2, 1.3, 1.1, 0.9)) +
    plot_annotation(
      title = sprintf("Figure 1 -- Head-mesoderm / bulge nuclei density (depth-integrated)%s", title_suffix),
      subtitle = sprintf(
        "Red contour = the data-driven head-mesoderm / bulge disc in each species.  Teal strip = context band used in panel C only (medaka landmark margin; zebrafish stable strip at theta %.0f° ± %d°).  0 min = ingression (medaka frame %d = %.1f min; zebrafish frame %d = %.1f min).  Density is always n_cells / area_um2 * 1000.",
        FIG1_ZEB_STRIP_CENTER, MARGIN_WIDTH_DEG,
        MEDAKA_INGRESSION_FRAME,
        FIG1_TIME_META[species == "Medaka", ingression_time_min],
        ZEB_INGRESSION_FRAME,
        FIG1_TIME_META[species == "Zebrafish", ingression_time_min]),
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(size = 9, color = "grey25")))
  save_pdf(fig, fname, w = 15, h = 21)

  # Per-region peak summary
  summary_dt <- d_ts[, {
    pk_raw <- peak_window_mean(time_min, density_per_1000um2)
    pk_align <- peak_window_mean(time_from_ingression_min,
                                 density_per_1000um2)
    dens_pre <- density_per_1000um2[
      time_from_ingression_min >= FIG1_BASELINE_FROM_INGR_MIN &
      time_from_ingression_min < 0
    ]
    if (!length(dens_pre)) dens_pre <- density_per_1000um2[time_from_ingression_min < 0]
    if (!length(dens_pre)) dens_pre <- head(density_per_1000um2, 5)
    .(pre_ingression_mean = round(mean(dens_pre, na.rm = TRUE), 2),
      peak_mean  = round(pk_raw$peak_mean, 2),
      peak_time  = round(pk_raw$peak_time),
      peak_time_from_ingression = round(pk_align$peak_time))
  }, by = .(species, region)]
  summary_dt[, peak_vs_pre_ingression := round(peak_mean / pre_ingression_mean, 2)]
  print(summary_dt)
  fwrite(d_ts, file.path(OUT_DIR,
                          sub("\\.pdf$", "_per_frame.csv", fname)))
  fwrite(ratio_ts, file.path(OUT_DIR,
                             sub("\\.pdf$", "_aligned_ratio.csv", fname)))
  fwrite(ratio_plot_dt, file.path(OUT_DIR,
                                  sub("\\.pdf$", "_aligned_ratio_binned.csv", fname)))
  list(summary = summary_dt, ratio = ratio_summary, density_qc = density_qc)
}

fig1_full <- build_fig1(
  "01_margin_density.pdf",
  time_caps = c(Medaka = Inf, Zebrafish = FIG1_ZEB_TCAP_MIN),
  title_suffix = sprintf("  (zebrafish t <= %d min)", FIG1_ZEB_TCAP_MIN))
fig1_t170 <- build_fig1(
  "01b_margin_density_t170.pdf",
  time_caps = c(Medaka = FIG1_ZEB_TCAP_MIN, Zebrafish = FIG1_ZEB_TCAP_MIN),
  title_suffix = sprintf("  (both species t <= %d min)", FIG1_ZEB_TCAP_MIN))

# =============================================================================
# FIGURE 2: MOVEMENT DIRECTION OVER TIME (both species)
# =============================================================================

banner("FIGURE 2 — movement direction over time")

TIME_BIN <- 30
vel_all <- rbind(
  vel_m[, .(species, time_min, v_epiboly, v_convergence, step_type,
            movement_type)],
  vel_z[, .(species, time_min, v_epiboly, v_convergence, step_type,
            movement_type)]
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

step_comp <- vel_all[!is.na(movement_type),
  .(n = .N), by = .(species, time_bin, movement_type)]
step_comp[, pct := 100 * n / sum(n), by = .(species, time_bin)]
step_comp[, movement_type := factor(movement_type, levels = STEP_TYPE_LEVELS)]

p2d <- ggplot(step_comp, aes(time_bin, pct, fill = movement_type)) +
  geom_area(alpha = 0.95, color = "white", linewidth = 0.2) +
  facet_wrap(~ species, ncol = 1) +
  scale_fill_manual(values = step_type_cols, drop = FALSE) +
  labs(title = "Movement type fraction over time (per-track classification)",
       subtitle = "Each step inherits its TRACK's movement type (same scheme as medaka_dynamics_05_spatial_flow).",
       x = "Time (min)", y = "% of steps", fill = NULL) +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 8, color = "grey25"))

fig2_caption <- paste0(
  "Step velocities (\u00b5m/min on sphere R):  ",
  "v_epiboly = R\u00b7d\u03b8/dt (>0 = vegetal-ward / epiboly direction)   |   ",
  "v_conv = R\u00b7sin\u03b8\u00b7d\u03c6/dt (along the equator; sign depends on flank)   |   ",
  "v_conv_dor = \u2212v_conv\u00b7sign(\u03c6\u2212\u03c6_bulge) (>0 = moving toward the dorsal axis = convergence)   |   ",
  "v_depth = d(spherical depth)/dt (>0 = moving inward = ingression).\n",
  "Per-TRACK classification (priority order):  ",
  sprintf("Ingression \u2014 track ENDS in ingression region AND total depth gain > P90 (medaka P90=%.1f \u00b5m within %.1f\u00b0 bulge disc; zebrafish P90=%.1f \u00b5m within margin strip \u00b1%g\u00b0).   ",
          trk_m$depth_thresh, BULGE_M$radius, trk_z$depth_thresh, INGR_STRIP_HALFW_Z),
  "Convergence \u2014 |mean v_conv_dor| > |mean v_epiboly| AND mean v_conv_dor > 0.   ",
  "Epiboly \u2014 mean v_epiboly > 0.   Animalward \u2014 otherwise.")

fig2 <- (p2a + p2b) / (p2c + p2d) +
  plot_annotation(title = "Movement direction over time",
                  caption = fig2_caption,
                  theme = theme(plot.title = element_text(face = "bold",
                                                          size = 14),
                                plot.caption = element_text(size = 8,
                                                             hjust = 0,
                                                             color = "grey25",
                                                             margin = margin(t = 10))))
save_pdf(fig2, "02_movement_direction.pdf", w = 16, h = 12)

fwrite(dir_summary, file.path(OUT_DIR, "direction_summary.csv"))
fwrite(step_comp, file.path(OUT_DIR, "step_type_composition.csv"))

# =============================================================================
# FIGURE 3: THICKNESS — both species in same plot
# Top:    thickness metric (IQR of SPHERICAL_DEPTH) per zone over time
# Bottom: depth distribution that the IQR summarises
#         (p10–p90 outer ribbon, p25–p75 inner ribbon = the IQR, median line)
# =============================================================================

banner("FIGURE 3 — thickness")

# Simplified zones: only three biologically meaningful tiers --
#   Animal cap    : tissue above the margin band (cells not yet involuted)
#   Margin        : whole margin strip (the entire involuting ring)
#   Head mesoderm : the data-driven deep-cell cluster (called "germ-ring
#                   bulge" in medaka and "shield" in zebrafish; these are
#                   biologically the same structure -- the site of cell
#                   internalisation at the dorsal margin).  Applied to BOTH
#                   species; the override replaces the local zone wherever
#                   a cell falls inside the data-driven disc.
ZONE_MARG_W_M  <- 5    # medaka  margin half-width (deg)
ZONE_MARG_W_Z  <- 8    # zebrafish margin half-width (deg)

ZONE_LEVELS <- c("Animal cap", "Margin", "Head mesoderm")

assign_zone <- function(theta, phi, margin_theta, marg_w, bulge) {
  zm_top <- margin_theta - marg_w
  z <- fcase(
    theta > zm_top,       "Margin",
    default               = "Animal cap")
  in_bulge <- !is.na(phi) &
    sqrt((theta - bulge$theta)^2 + (phi - bulge$phi)^2) < bulge$radius
  z[in_bulge] <- "Head mesoderm"   # overrides any local zone
  factor(z, levels = ZONE_LEVELS)
}

sp_m[, zone := assign_zone(THETA_DEG, PHI_DEG, MARGIN_M,
                            ZONE_MARG_W_M, BULGE_M)]
sp_z[, zone := assign_zone(THETA_DEG, PHI_DEG, MARGIN_Z,
                            ZONE_MARG_W_Z, BULGE_Z)]

# Per-cell depth distribution in 10-min bins per zone, per species.
# Built ONCE on the full data; the t<=170 variant of Fig 3 just filters this.
depth_dist <- rbind(
  sp_m[!is.na(zone) & is.finite(SPHERICAL_DEPTH),
       .(time_bin = floor(time_min / 10) * 10,
         zone, depth_um = SPHERICAL_DEPTH, species = "Medaka")],
  sp_z[!is.na(zone) & is.finite(SPHERICAL_DEPTH),
       .(time_bin = floor(time_min / 10) * 10,
         zone, depth_um = SPHERICAL_DEPTH, species = "Zebrafish")]
)
depth_dist[, species := factor(species, levels = c("Medaka", "Zebrafish"))]

BIN_DEG <- 5
make_window <- function(sp_dt, sp_label, t_cap = Inf) {
  d <- sp_dt[time_min <= t_cap]
  tmax <- max(d$time_min, na.rm = TRUE)
  wlen <- 30                      # window width in minutes
  mid_t <- tmax / 2
  rbind(
    d[time_min <= min(wlen, tmax),
      .(species = sp_label,
        window  = sprintf("Early\n(0\u2013%.0f min)", min(wlen, tmax)),
        window_order = 1L,
        THETA_DEG, PHI_DEG, SPHERICAL_DEPTH)],
    d[time_min >= mid_t - wlen / 2 & time_min <= mid_t + wlen / 2,
      .(species = sp_label,
        window  = sprintf("Mid\n(%.0f\u2013%.0f min)",
                          mid_t - wlen / 2, mid_t + wlen / 2),
        window_order = 2L,
        THETA_DEG, PHI_DEG, SPHERICAL_DEPTH)],
    d[time_min >= max(0, tmax - wlen),
      .(species = sp_label,
        window  = sprintf("Late\n(%.0f\u2013%.0f min)", max(0, tmax - wlen), tmax),
        window_order = 3L,
        THETA_DEG, PHI_DEG, SPHERICAL_DEPTH)]
  )
}

# -----------------------------------------------------------------------------
# Fig 3 builder.  Two versions are produced:
#   1) full data  -> 03_thickness.pdf
#   2) t <= 170 min for both species -> 03b_thickness_t170.pdf
#
# Panel A metric is the p10-p90 spherical depth range of the cells in each
# zone at each time bin (the width of the lighter ribbon in panel B).
# Panel B shows that p10-p90 ribbon and the IQR (p25-p75) ribbon + median.
# -----------------------------------------------------------------------------
build_fig3 <- function(t_cap, fname, title_suffix = "") {
  dd <- depth_dist[time_bin <= t_cap]
  ds <- dd[, .(p10 = quantile(depth_um, 0.10),
               p25 = quantile(depth_um, 0.25),
               p50 = quantile(depth_um, 0.50),
               p75 = quantile(depth_um, 0.75),
               p90 = quantile(depth_um, 0.90),
               full_range = quantile(depth_um, 0.90) - quantile(depth_um, 0.10),
               n_cells = .N),
           by = .(species, zone, time_bin)][n_cells >= 20]

  p3a <- ggplot(ds, aes(time_bin, full_range, color = species)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.4, alpha = 0.7) +
    facet_wrap(~ zone, nrow = 1, scales = "free_x") +
    scale_color_manual(values = species_colors) +
    labs(title = "Tissue thickness \u2014 p10\u2013p90 spherical depth range per zone",
         x = "Time (min)", y = "thickness = p90 \u2212 p10 depth (\u00b5m)",
         color = NULL) +
    theme_pub()

  p3b <- ggplot(ds, aes(x = time_bin)) +
    geom_ribbon(aes(ymin = p10,  ymax = p90,  fill = species), alpha = 0.20) +
    geom_ribbon(aes(ymin = p25,  ymax = p75,  fill = species), alpha = 0.40) +
    geom_line(aes(y = p50, color = species), linewidth = 0.9) +
    facet_grid(species ~ zone, scales = "free") +
    scale_fill_manual(values = species_colors, guide = "none") +
    scale_color_manual(values = species_colors, guide = "none") +
    labs(title = "Depth distribution of cells per zone (the data underlying the thickness above)",
         subtitle = "Lighter band = p10\u2013p90 (= thickness above).  Darker band = IQR.  Line = median.  y-axis: higher = deeper.",
         x = "Time (min)", y = expression("spherical depth ("*mu*"m; higher = deeper)")) +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 8, color = "grey25"))

  heat_raw <- rbind(make_window(sp_m, "Medaka",    t_cap),
                    make_window(sp_z, "Zebrafish", t_cap))
  heat_raw[, phi_bin   := floor(PHI_DEG   / BIN_DEG) * BIN_DEG + BIN_DEG / 2]
  heat_raw[, theta_bin := floor(THETA_DEG / BIN_DEG) * BIN_DEG + BIN_DEG / 2]
  heat <- heat_raw[is.finite(SPHERICAL_DEPTH),
    .(deep_p90 = quantile(SPHERICAL_DEPTH, 0.90, na.rm = TRUE), n = .N),
    by = .(species, window, window_order, phi_bin, theta_bin)][n >= 3]
  heat[, panel_col := factor(c("Early", "Mid", "Late")[window_order],
                              levels = c("Early", "Mid", "Late"))]
  heat[, species := factor(species, levels = c("Medaka", "Zebrafish"))]
  tmax_m <- min(t_cap, max(sp_m$time_min, na.rm = TRUE))
  tmax_z <- min(t_cap, max(sp_z$time_min, na.rm = TRUE))

  heat_band <- data.table(
    species = factor(c("Medaka", "Zebrafish"), levels = c("Medaka", "Zebrafish")),
    ymin    = c(MARGIN_M - MARGIN_WIDTH_DEG, MARGIN_Z - MARGIN_WIDTH_DEG),
    ymax    = c(MARGIN_M + MARGIN_WIDTH_DEG, MARGIN_Z + MARGIN_WIDTH_DEG))

  p3c <- ggplot(heat, aes(phi_bin, theta_bin, fill = deep_p90)) +
    geom_tile() +
    geom_rect(data = heat_band, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              fill = NA, color = "white", linewidth = 0.4, linetype = "dashed") +
    geom_polygon(data = bulge_poly, inherit.aes = FALSE,
                 aes(phi, theta), fill = NA,
                 color = "white", linewidth = 0.5) +
    geom_text(data = unique(heat[, .(species, panel_col, window)]),
              aes(x = -Inf, y = -Inf, label = window),
              inherit.aes = FALSE, hjust = -0.05, vjust = -0.5,
              size = 2.6, color = "grey80", lineheight = 0.85) +
    facet_grid(species ~ panel_col, scales = "free") +
    scale_fill_viridis_c(option = "inferno", name = "P90 depth (\u00b5m)",
                         limits = c(0, NA), oob = scales::squish) +
    scale_y_reverse() +
    labs(title = "Where the tissue is thickest on the embryo surface (temporal evolution)",
         subtitle = sprintf("Tile = %g\u00b0 \u00d7 %g\u00b0 (\u03c6 \u00d7 \u03b8), coloured by P90 spherical depth.  Dashed = margin band.  Solid ellipse = head-mesoderm disc.  Medaka 0\u2013%.0f min, zebrafish 0\u2013%.0f min.",
                            BIN_DEG, BIN_DEG, tmax_m, tmax_z),
         x = expression(varphi*" (deg)"), y = expression(theta*" (deg, animal pole up)")) +
    theme_pub() +
    theme(plot.subtitle = element_text(size = 9, color = "grey25"))

  fig <- (p3a / p3b / p3c) +
    plot_layout(heights = c(1, 1.6, 1.6)) +
    plot_annotation(
      title = sprintf("Tissue thickness over time \u2014 Medaka vs Zebrafish%s",
                      title_suffix),
      subtitle = sprintf("Thickness = p90 \u2212 p10 spherical depth per zone.  Zones: Animal cap = \u03b8 \u2264 margin-W\u00b0 ; Margin = below that ; Head mesoderm = bulge disc (overrides). W = %g\u00b0 medaka, %g\u00b0 zebrafish.",
                          ZONE_MARG_W_M, ZONE_MARG_W_Z),
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(size = 9, color = "grey25")))
  save_pdf(fig, fname, w = 17, h = 18)
  ds
}

depth_summary       <- build_fig3(Inf, "03_thickness.pdf")
depth_summary_t170  <- build_fig3(170, "03b_thickness_t170.pdf",
                                   title_suffix = "  (t \u2264 170 min)")

fwrite(depth_summary,      file.path(OUT_DIR, "thickness_depth_distribution.csv"))
fwrite(depth_summary_t170, file.path(OUT_DIR, "thickness_depth_distribution_t170.csv"))

# Quick numeric check -- does the margin layer actually thicken?
zone_report <- list(
  Medaka    = c("Margin", "Head mesoderm"),
  Zebrafish = c("Margin", "Head mesoderm"))
for (sp in c("Medaka", "Zebrafish")) {
  cat(sprintf("  --- %s ---\n", sp))
  for (zn in zone_report[[sp]]) {
    dd <- depth_summary[zone == zn & species == sp][order(time_bin)]
    if (nrow(dd) >= 3) {
      t0   <- dd$full_range[1]
      tmax <- max(dd$full_range)
      tend <- tail(dd$full_range, 1)
      tmax_t <- dd$time_bin[which.max(dd$full_range)]
      cat(sprintf("    %-26s full-range: t=0 %.1f \u00b5m | peak %.1f \u00b5m at t=%g min | end %.1f \u00b5m (\u0394peak=%+.1f%%, \u0394end=%+.1f%%)\n",
                  zn, t0, tmax, tmax_t, tend,
                  100 * (tmax - t0) / t0, 100 * (tend - t0) / t0))
    }
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

# -----------------------------------------------------------------------------
# 04d: Ingression-only vector map (validation of the per-TRACK classifier)
# -----------------------------------------------------------------------------
# Uses the SAME track-level movement_type produced by classify_tracks() that
# Figure 2 uses (i.e. nothing here re-classifies anything). For each track
# tagged "Ingression" we draw a SINGLE net-displacement arrow from its start
# to its end position, so the per-step jitter does not turn the panel into
# a fur-ball.  Two views, faceted by the HOUR IN WHICH THE TRACK STARTS:
#
#   * Top view (XY)  — shows WHERE on the embryo surface ingression happens
#                      (should cluster inside the medaka bulge ellipse and
#                      along the zebrafish margin strip).
#   * Side view (YZ) — shows the cells actually diving INWARD (this is the
#                      signal the classifier keys on: depth gain).  Arrows
#                      pointing toward the embryo interior validate that
#                      "Ingression" tracks really move into the deep tissue.
#
# Arrow colour = depth change (last - first SPHERICAL_DEPTH, in µm).
# -----------------------------------------------------------------------------
banner("FIGURE 4d — ingression-only net-displacement vectors (classifier validation)")

ingression_track_endpoints <- function(vel, sp_label) {
  setkey(vel, TRACK_ID, FRAME)
  ing <- vel[movement_type == "Ingression",
             .(start_frame = FRAME[1],
               end_frame   = FRAME[.N],
               start_x = POSITION_X[1], end_x = POSITION_X[.N],
               start_y = POSITION_Y[1], end_y = POSITION_Y[.N],
               start_z = POSITION_Z[1], end_z = POSITION_Z[.N],
               start_depth = SPHERICAL_DEPTH[1],
               end_depth   = SPHERICAL_DEPTH[.N],
               depth_change = SPHERICAL_DEPTH[.N] - SPHERICAL_DEPTH[1],
               n_steps     = .N),
             by = TRACK_ID]
  ing[, start_hour := floor(start_frame * (if (sp_label == "Medaka") MEDAKA_FI
                                            else ZEB_FI) / 3600)]
  ing[, species := sp_label]
  ing
}
ing_ep_m <- ingression_track_endpoints(vel_m, "Medaka")
ing_ep_z <- ingression_track_endpoints(vel_z, "Zebrafish")

cat(sprintf("  Ingression tracks: Medaka n=%d (median depth gain %.1f µm) | Zebrafish n=%d (median depth gain %.1f µm)\n",
            nrow(ing_ep_m), median(ing_ep_m$depth_change),
            nrow(ing_ep_z), median(ing_ep_z$depth_change)))

# Per-hour facet labels
ing_hour_labels <- function(ep) {
  ep[, .(n_tracks = .N,
         med_depth_gain = round(median(depth_change), 1)),
     by = start_hour][order(start_hour)][
       , hour_lab := sprintf("h%02d  (n=%d, Δdepth med %+.1f µm)",
                             start_hour, n_tracks, med_depth_gain)][]
}
ihl_m <- ing_hour_labels(ing_ep_m)
ihl_z <- ing_hour_labels(ing_ep_z)

plot_ingression_endpoints <- function(ep, ihl, sp_label) {
  if (nrow(ep) == 0) return(NULL)
  ep <- merge(ep, ihl[, .(start_hour, hour_lab)], by = "start_hour")
  # Colour limits clipped to robust range
  dc_lim <- range(ep$depth_change, na.rm = TRUE)
  dc_lim[1] <- max(dc_lim[1], quantile(ep$depth_change, 0.02, na.rm = TRUE))
  dc_lim[2] <- min(dc_lim[2], quantile(ep$depth_change, 0.98, na.rm = TRUE))

  top <- ggplot(ep) +
    geom_segment(aes(x = start_x, y = start_y,
                     xend = end_x, yend = end_y,
                     color = depth_change),
                 arrow = arrow(length = unit(0.10, "cm"), type = "closed"),
                 linewidth = 0.5, alpha = 0.9) +
    geom_point(aes(start_x, start_y), size = 0.4,
               color = "grey30", alpha = 0.6) +
    facet_wrap(~ hour_lab, ncol = 4) +
    scale_color_viridis_c(option = "magma", direction = -1,
                          limits = dc_lim, oob = scales::squish,
                          name = "Δ depth (µm)\nfirst→last") +
    scale_y_reverse() + coord_fixed() +
    labs(title = sprintf("%s — INGRESSION tracks: TOP view (XY)  — WHERE on the surface", sp_label),
         subtitle = "one arrow per track = start→end XY position | dot = start | colour = depth gained over track lifetime",
         x = "X (µm)", y = "Y (µm)  [AP top, VP bottom]") +
    theme_pub(10) + theme(strip.text = element_text(size = 8))

  side <- ggplot(ep) +
    geom_segment(aes(x = start_y, y = start_z,
                     xend = end_y, yend = end_z,
                     color = depth_change),
                 arrow = arrow(length = unit(0.10, "cm"), type = "closed"),
                 linewidth = 0.5, alpha = 0.9) +
    geom_point(aes(start_y, start_z), size = 0.4,
               color = "grey30", alpha = 0.6) +
    facet_wrap(~ hour_lab, ncol = 4) +
    scale_color_viridis_c(option = "magma", direction = -1,
                          limits = dc_lim, oob = scales::squish,
                          name = "Δ depth (µm)\nfirst→last") +
    coord_fixed() +
    labs(title = sprintf("%s — INGRESSION tracks: SIDE view (YZ)  — diving INWARD", sp_label),
         subtitle = "arrows should point toward the embryo interior (Z away from the outer surface) = the signal the classifier keys on",
         x = "Y (µm)", y = "Z (µm)") +
    theme_pub(10) + theme(strip.text = element_text(size = 8))

  list(top = top, side = side)
}

# Per-track endpoint CSV (useful for downstream inspection)
fwrite(rbind(ing_ep_m, ing_ep_z, fill = TRUE),
       file.path(OUT_DIR, "ingression_track_endpoints.csv"))

ing_pl_m <- plot_ingression_endpoints(ing_ep_m, ihl_m, "Medaka")
ing_pl_z <- plot_ingression_endpoints(ing_ep_z, ihl_z, "Zebrafish")

if (!is.null(ing_pl_m)) {
  n_h <- nrow(ihl_m); h_panel <- max(6, 4 * ceiling(n_h / 4))
  save_pdf(ing_pl_m$top / ing_pl_m$side + plot_layout(heights = c(1, 1)),
           "04d_ingression_vectors_medaka.pdf", w = 18, h = 2 * h_panel)
}
if (!is.null(ing_pl_z)) {
  n_h <- nrow(ihl_z); h_panel <- max(6, 4 * ceiling(n_h / 4))
  save_pdf(ing_pl_z$top / ing_pl_z$side + plot_layout(heights = c(1, 1)),
           "04e_ingression_vectors_zebrafish.pdf", w = 18, h = 2 * h_panel)
}

fwrite(rbind(ihl_m[, species := "Medaka"], ihl_z[, species := "Zebrafish"],
             fill = TRUE),
       file.path(OUT_DIR, "ingression_tracks_per_hour.csv"))

# =============================================================================
# FIGURE 5: TRACK COMPARISON (FAIR SAMPLING)
# =============================================================================
# Filter: 20 <= duration <= 60 min, full-track net displacement >= 20 um.
# To compare two species sampled at different rates (medaka 30 s, zebrafish
# 120 s), we ALWAYS use a common 120-s, non-overlapping step interval here.
# Track duration is measured on the full track (sampling-invariant); path
# length, turning, and MSD all use the 120-s steps so neither species
# accumulates extra jitter.  The QC output is a lightweight HTML animation of a
# fair-sampled subset of these tracks over real time.
# =============================================================================

banner("FIGURE 5 — track comparison (fair 120-s sampling)")

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

matched_endpoint <- ep_all[duration_min >= FIG5_MIN_DURATION_MIN &
    duration_min <= FIG5_MAX_DURATION_MIN &
    net_disp_um >= FIG5_MIN_NET_DISP_UM]
ids_m <- matched_endpoint[species == "Medaka",    TRACK_ID]
ids_z <- matched_endpoint[species == "Zebrafish", TRACK_ID]
cat(sprintf("  Duration-filtered tracks (%d-%d min, >= %d um endpoint displacement): Medaka n=%d, Zebrafish n=%d\n",
       FIG5_MIN_DURATION_MIN, FIG5_MAX_DURATION_MIN,
  FIG5_MIN_NET_DISP_UM,
       length(ids_m), length(ids_z)))
fwrite(matched_endpoint,
  file.path(OUT_DIR, "duration_filtered_track_endpoints.csv"))

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
fwrite(matched, file.path(OUT_DIR, "duration_filtered_tracks_20to60min.csv"))

metric_ids_m <- matched[species == "Medaka", TRACK_ID]
metric_ids_z <- matched[species == "Zebrafish", TRACK_ID]

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

track_surface_qc <- function(sp, ids, sp_label) {
  if (!length(ids)) return(data.table())
  d <- copy(sp[TRACK_ID %in% ids])
  setorder(d, TRACK_ID, FRAME)
  d[, .(
    start_phi = PHI_DEG[1],
    start_theta = THETA_DEG[1],
    end_phi = PHI_DEG[.N],
    end_theta = THETA_DEG[.N]
  ), by = TRACK_ID][, species := sp_label][]
}

track_path_qc <- function(sub, ids, sp_label,
                          fi_sec,
                          max_tracks = FIG5_QC_SAMPLE_TRACKS,
                          seed = FIG5_QC_SEED) {
  if (!length(ids)) {
    return(list(paths = data.table(), starts = data.table(), ends = data.table()))
  }
  keep_ids <- ids
  if (length(keep_ids) > max_tracks) {
    set.seed(seed)
    keep_ids <- sample(keep_ids, max_tracks)
  }
  d <- copy(sub[TRACK_ID %in% keep_ids,
                .(TRACK_ID, FRAME, POSITION_X, POSITION_Y)])
  setorder(d, TRACK_ID, FRAME)
  d[, `:=`(species = sp_label,
           time_min = FRAME * fi_sec / 60)]
  se <- d[, .(
    start_x = POSITION_X[1],
    start_y = POSITION_Y[1],
    end_x   = POSITION_X[.N],
    end_y   = POSITION_Y[.N]
  ), by = TRACK_ID]
  se[, species := sp_label]
  list(
    paths = d,
    starts = se[, .(species, TRACK_ID, x = start_x, y = start_y)],
    ends = se[, .(species, TRACK_ID, x = end_x, y = end_y)]
  )
}

build_track_qc_payload <- function(path_dt) {
  d <- copy(path_dt)
  d[, species := as.character(species)]
  setorder(d, species, TRACK_ID, time_min)
  time_values <- sort(unique(d$time_min))
  species_payload <- setNames(vector("list", 2L), c("Medaka", "Zebrafish"))
  for (sp_label in c("Medaka", "Zebrafish")) {
    sp <- d[species == sp_label]
    if (!nrow(sp)) {
      species_payload[[sp_label]] <- list(
        label = sp_label,
        color = unname(species_colors[sp_label]),
        n_tracks = 0L,
        xmin = 0,
        xmax = 1,
        ymin = 0,
        ymax = 1,
        tracks = list()
      )
      next
    }
    track_list <- sp[, .(
      track = list(list(
        id = as.character(TRACK_ID[1]),
        time_min = time_min,
        x = POSITION_X,
        y = POSITION_Y
      ))
    ), by = TRACK_ID]$track
    x_rng <- range(sp$POSITION_X, na.rm = TRUE)
    y_rng <- range(sp$POSITION_Y, na.rm = TRUE)
    species_payload[[sp_label]] <- list(
      label = sp_label,
      color = unname(species_colors[sp_label]),
      n_tracks = length(track_list),
      xmin = x_rng[1],
      xmax = x_rng[2],
      ymin = y_rng[1],
      ymax = y_rng[2],
      tracks = track_list
    )
  }
  list(
    times = time_values,
    trail_steps = FIG5_QC_TRAIL_STEPS,
    species = species_payload
  )
}

build_track_qc_html <- function(payload) {
  payload_json <- jsonlite::toJSON(payload, auto_unbox = TRUE,
                                   digits = 7, null = "null")
  subtitle <- sprintf(
    "Tracks: %d-%d min and >= %d um full-track endpoint displacement. Random sample <= %d tracks per species. Each dot is the current position and each tail shows the last %d fair-sampled steps.",
    FIG5_MIN_DURATION_MIN,
    FIG5_MAX_DURATION_MIN,
    FIG5_MIN_NET_DISP_UM,
    FIG5_QC_SAMPLE_TRACKS,
    FIG5_QC_TRAIL_STEPS
  )
  c(
    "<!DOCTYPE html>",
    "<html lang=\"en\">",
    "<head>",
    "  <meta charset=\"utf-8\">",
    "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "  <title>Duration-filtered track QC viewer</title>",
    "  <style>",
    "    :root { color-scheme: light; font-family: Helvetica, Arial, sans-serif; }",
    "    body { margin: 0; background: #f6f3eb; color: #1f1f1f; }",
    "    .page { max-width: 1500px; margin: 0 auto; padding: 24px; }",
    "    h1 { margin: 0 0 8px; font-size: 28px; }",
    "    .subtitle { margin: 0 0 20px; max-width: 1050px; line-height: 1.45; color: #4b4b4b; }",
    "    .controls { display: flex; gap: 12px; align-items: center; flex-wrap: wrap; margin-bottom: 18px; }",
    "    button { border: 0; border-radius: 999px; background: #1f1f1f; color: #ffffff; padding: 10px 18px; font-size: 15px; cursor: pointer; }",
    "    input[type=range] { flex: 1 1 380px; accent-color: #1f1f1f; }",
    "    .time-label { min-width: 108px; font-variant-numeric: tabular-nums; }",
    "    .grid { display: grid; grid-template-columns: repeat(2, minmax(320px, 1fr)); gap: 18px; }",
    "    .panel { background: #fffdf8; border: 1px solid #ddd4c6; border-radius: 18px; padding: 16px; box-shadow: 0 12px 28px rgba(43, 32, 17, 0.08); }",
    "    .panel h2 { margin: 0 0 8px; font-size: 18px; }",
    "    .meta { margin: 0 0 10px; color: #6a6256; font-size: 14px; }",
    "    canvas { width: 100%; height: auto; display: block; background: #fffdf8; border-radius: 12px; }",
    "    .footer { margin-top: 12px; color: #6a6256; font-size: 13px; }",
    "    @media (max-width: 960px) { .grid { grid-template-columns: 1fr; } }",
    "  </style>",
    "</head>",
    "<body>",
    "  <div class=\"page\">",
    "    <h1>Figure 5 QC viewer</h1>",
    paste0("    <p class=\"subtitle\">", subtitle, "</p>"),
    "    <div class=\"controls\">",
    "      <button id=\"playButton\" type=\"button\">Play</button>",
    "      <input id=\"timeSlider\" type=\"range\" min=\"0\" max=\"0\" step=\"1\" value=\"0\">",
    "      <div id=\"timeLabel\" class=\"time-label\">Time: 0 min</div>",
    "    </div>",
    "    <div class=\"grid\">",
    "      <section class=\"panel\">",
    "        <h2>Medaka</h2>",
    "        <p id=\"metaMedaka\" class=\"meta\"></p>",
    "        <canvas id=\"canvasMedaka\" width=\"900\" height=\"900\"></canvas>",
    "      </section>",
    "      <section class=\"panel\">",
    "        <h2>Zebrafish</h2>",
    "        <p id=\"metaZebrafish\" class=\"meta\"></p>",
    "        <canvas id=\"canvasZebrafish\" width=\"900\" height=\"900\"></canvas>",
    "      </section>",
    "    </div>",
    "    <p class=\"footer\">The viewer uses the same 120 s fair-sampled tracks used for straightness, turning angle, and MSD.</p>",
    "  </div>",
    "  <script>",
    paste0("    const qcData = ", payload_json, ";"),
    "    const slider = document.getElementById('timeSlider');",
    "    const playButton = document.getElementById('playButton');",
    "    const timeLabel = document.getElementById('timeLabel');",
    "    const panelNames = ['Medaka', 'Zebrafish'];",
    "    const canvasIds = { Medaka: 'canvasMedaka', Zebrafish: 'canvasZebrafish' };",
    "    const metaIds = { Medaka: 'metaMedaka', Zebrafish: 'metaZebrafish' };",
    "    const padding = { left: 58, right: 18, top: 18, bottom: 52 };",
    "    let frameIndex = 0;",
    "    let timerId = null;",
    "    const states = {};",
    "",
    "    function init() {",
    "      panelNames.forEach((name) => {",
    "        const panel = qcData.species[name];",
    "        const canvas = document.getElementById(canvasIds[name]);",
    "        states[name] = { panel, canvas, ctx: canvas.getContext('2d') };",
    "        document.getElementById(metaIds[name]).textContent = `${panel.n_tracks} sampled tracks`;",
    "      });",
    "      slider.max = Math.max(qcData.times.length - 1, 0);",
    "      slider.addEventListener('input', (event) => {",
    "        frameIndex = Number(event.target.value);",
    "        render();",
    "      });",
    "      playButton.addEventListener('click', () => {",
    "        if (timerId === null) {",
    "          if (frameIndex >= qcData.times.length - 1) frameIndex = 0;",
    "          play();",
    "        } else {",
    "          pause();",
    "        }",
    "      });",
    "      render();",
    "    }",
    "",
    "    function play() {",
    "      playButton.textContent = 'Pause';",
    "      timerId = window.setInterval(() => {",
    "        if (frameIndex >= qcData.times.length - 1) {",
    "          pause();",
    "          return;",
    "        }",
    "        frameIndex += 1;",
    "        render();",
    "      }, 360);",
    "    }",
    "",
    "    function pause() {",
    "      if (timerId !== null) {",
    "        window.clearInterval(timerId);",
    "        timerId = null;",
    "      }",
    "      playButton.textContent = 'Play';",
    "    }",
    "",
    "    function xToCanvas(state, x) {",
    "      const innerWidth = state.canvas.width - padding.left - padding.right;",
    "      const span = Math.max(state.panel.xmax - state.panel.xmin, 1e-6);",
    "      return padding.left + ((x - state.panel.xmin) / span) * innerWidth;",
    "    }",
    "",
    "    function yToCanvas(state, y) {",
    "      const innerHeight = state.canvas.height - padding.top - padding.bottom;",
    "      const span = Math.max(state.panel.ymax - state.panel.ymin, 1e-6);",
    "      return padding.top + ((y - state.panel.ymin) / span) * innerHeight;",
    "    }",
    "",
    "    function lastVisibleIndex(track, currentTime) {",
    "      for (let i = track.time_min.length - 1; i >= 0; i -= 1) {",
    "        if (track.time_min[i] <= currentTime + 1e-9) return i;",
    "      }",
    "      return -1;",
    "    }",
    "",
    "    function drawAxes(state) {",
    "      const ctx = state.ctx;",
    "      const width = state.canvas.width;",
    "      const height = state.canvas.height;",
    "      ctx.fillStyle = '#fffdf8';",
    "      ctx.fillRect(0, 0, width, height);",
    "      ctx.strokeStyle = '#ddd4c6';",
    "      ctx.lineWidth = 1.5;",
    "      ctx.strokeRect(padding.left, padding.top, width - padding.left - padding.right, height - padding.top - padding.bottom);",
    "      ctx.fillStyle = '#6a6256';",
    "      ctx.font = '24px Helvetica, Arial, sans-serif';",
    "      ctx.fillText('X (um)', width / 2 - 34, height - 12);",
    "      ctx.save();",
    "      ctx.translate(18, height / 2 + 34);",
    "      ctx.rotate(-Math.PI / 2);",
    "      ctx.fillText('Y (um)', 0, 0);",
    "      ctx.restore();",
    "    }",
    "",
    "    function drawTracks(state, currentTime) {",
    "      const ctx = state.ctx;",
    "      const color = state.panel.color;",
    "      state.panel.tracks.forEach((track) => {",
    "        const idx = lastVisibleIndex(track, currentTime);",
    "        if (idx < 0) return;",
    "        const firstIdx = Math.max(0, idx - qcData.trail_steps + 1);",
    "        ctx.beginPath();",
    "        for (let i = firstIdx; i <= idx; i += 1) {",
    "          const px = xToCanvas(state, track.x[i]);",
    "          const py = yToCanvas(state, track.y[i]);",
    "          if (i === firstIdx) ctx.moveTo(px, py);",
    "          else ctx.lineTo(px, py);",
    "        }",
    "        ctx.strokeStyle = color;",
    "        ctx.lineWidth = 1.25;",
    "        ctx.globalAlpha = 0.28;",
    "        ctx.stroke();",
    "        const cx = xToCanvas(state, track.x[idx]);",
    "        const cy = yToCanvas(state, track.y[idx]);",
    "        ctx.globalAlpha = 0.95;",
    "        ctx.fillStyle = color;",
    "        ctx.beginPath();",
    "        ctx.arc(cx, cy, 4.0, 0, Math.PI * 2);",
    "        ctx.fill();",
    "        ctx.strokeStyle = '#111111';",
    "        ctx.lineWidth = 0.8;",
    "        ctx.stroke();",
    "      });",
    "      ctx.globalAlpha = 1;",
    "    }",
    "",
    "    function render() {",
    "      if (!qcData.times.length) return;",
    "      const currentTime = qcData.times[frameIndex];",
    "      timeLabel.textContent = `Time: ${currentTime.toFixed(0)} min`;",
    "      slider.value = frameIndex;",
    "      panelNames.forEach((name) => {",
    "        drawAxes(states[name]);",
    "        drawTracks(states[name], currentTime);",
    "      });",
    "    }",
    "",
    "    init();",
    "  </script>",
    "</body>",
    "</html>"
  )
}

surface_qc <- rbind(
  track_surface_qc(sp_m, metric_ids_m, "Medaka"),
  track_surface_qc(sp_z, metric_ids_z, "Zebrafish")
)
surface_qc[, species := factor(species, levels = c("Medaka", "Zebrafish"))]

path_qc_m <- track_path_qc(sub_m, metric_ids_m, "Medaka", MEDAKA_FI)
path_qc_z <- track_path_qc(sub_z, metric_ids_z, "Zebrafish", ZEB_FI)
paths_qc <- rbindlist(list(path_qc_m$paths, path_qc_z$paths), fill = TRUE)
paths_qc[, species := factor(species, levels = c("Medaka", "Zebrafish"))]

fwrite(surface_qc, file.path(OUT_DIR, "duration_filtered_track_starts.csv"))
fwrite(paths_qc, file.path(OUT_DIR, "duration_filtered_track_sample_paths.csv"))
old_qc_pdf <- file.path(OUT_DIR, "05a_duration_filtered_track_qc.pdf")
if (file.exists(old_qc_pdf)) unlink(old_qc_pdf)
qc_payload <- build_track_qc_payload(paths_qc)
qc_html <- build_track_qc_html(qc_payload)
save_html_file(qc_html, "05a_duration_filtered_track_qc.html")

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
msd_m <- msd_matched_fair(sp_m, metric_ids_m, MEDAKA_FI, 4L, max_k)[, species := "Medaka"]
msd_z <- msd_matched_fair(sp_z, metric_ids_z, ZEB_FI,    1L, max_k)[, species := "Zebrafish"]
msd <- rbind(msd_m, msd_z)
fwrite(msd, file.path(OUT_DIR, "matched_msd.csv"))
fwrite(msd, file.path(OUT_DIR, "duration_filtered_msd_20to60min.csv"))

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
      subtitle = sprintf("20-60 min duration, >=20 um endpoint displacement, with valid fair-sampled metrics (Medaka n=%d, Zebrafish n=%d) | common 120-s lags",
             length(metric_ids_m), length(metric_ids_z)),
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
    title = "Track comparison (fair 120-s sampling)",
    subtitle = sprintf("Tracks: %d-%d min full-track duration and >= %d um endpoint displacement. Path length, turning, and MSD all use common 120-s steps.",
                       FIG5_MIN_DURATION_MIN, FIG5_MAX_DURATION_MIN,
                       FIG5_MIN_NET_DISP_UM),
    theme = theme(plot.title = element_text(face = "bold", size = 14)))
save_pdf(fig5, "05_matched_tracks.pdf", w = 13, h = 11)

banner("DONE")
cat(sprintf("  outputs in %s/\n", OUT_DIR))
