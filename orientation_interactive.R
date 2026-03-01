# =============================================================================
# Enhanced Interactive Embryo Orientation, Sphere Fitting & Margin Selection
# =============================================================================
#
# PURPOSE:
#   Step 1 of the comparative medaka/zebrafish analysis pipeline.
#   Provides interactive 3D visualisation and annotation of light-sheet
#   microscopy nuclear-tracking data (TrackMate exports).
#
# FEATURES:
#   1. Fast loading of large TrackMate CSVs (millions of spots via data.table)
#   2. Temporal navigation – scrub through developmental time with track tails
#   3. Embryo orientation (animal pole → +Y, dorsal → +X)
#   4. Robust sphere fitting (works with hemispheres / partial embryo data)
#   5. Spherical coordinate system (depth from surface, latitude, longitude)
#   6. Margin altitude selection for ingression analysis & cell flagging
#   7. Export enriched CSVs with spherical coords + margin flags
#
# OUTPUT (downloaded as ZIP):
#   - oriented_spots.csv   — all spots with spherical coords & IN_MARGIN flag
#   - oriented_tracks.csv  — tracks with reoriented coordinates
#   - sphere_params.csv    — fitted sphere centre, radius, margin bounds
#
# PIPELINE:
#   → Step 1: THIS SCRIPT
#     Step 2: trackmate_filter_and_validate.R
#     Step 3: trackmate_analysis.R
#
# =============================================================================

library(dplyr)
library(readr)
library(data.table)
library(plotly)
library(shiny)
library(pracma)

options(shiny.maxRequestSize = 2000 * 1024^2)  # 2 GB upload limit

# #############################################################################
#                        CORE FUNCTIONS — Data Loading
# #############################################################################

#' Fast load of TrackMate CSVs.
#' Uses data.table::fread for speed; auto-detects and skips TM metadata rows.
load_trackmate_data <- function(spots_file, tracks_file = NULL) {
  spots_raw <- as_tibble(data.table::fread(spots_file, header = TRUE,
                                           showProgress = FALSE))
  # Detect TrackMate 3-row metadata
  if (nrow(spots_raw) > 3 && "POSITION_X" %in% names(spots_raw)) {
    if (is.na(suppressWarnings(as.numeric(spots_raw$POSITION_X[1])))) {
      spots_raw <- spots_raw[-c(1:3), ]
    }
  }

  spots <- spots_raw %>%
    mutate(across(any_of(c("ID", "TRACK_ID", "POSITION_X", "POSITION_Y",
                           "POSITION_Z", "FRAME", "QUALITY")),
                  ~suppressWarnings(as.numeric(.)))) %>%
    filter(!is.na(TRACK_ID))

  tracks <- NULL
  if (!is.null(tracks_file)) {
    tracks_raw <- as_tibble(data.table::fread(tracks_file, header = TRUE,
                                              showProgress = FALSE))
    if (nrow(tracks_raw) > 3) {
      tc <- intersect(names(tracks_raw),
                      c("TRACK_DURATION", "NUMBER_SPOTS", "TRACK_MEAN_SPEED"))
      if (length(tc) > 0 &&
          is.na(suppressWarnings(as.numeric(tracks_raw[[tc[1]]][1])))) {
        tracks_raw <- tracks_raw[-c(1:3), ]
      }
    }
    tracks <- tracks_raw %>%
      mutate(across(-any_of("LABEL"), ~suppressWarnings(as.numeric(.))))
  }

  list(spots = spots, tracks = tracks)
}

# #############################################################################
#                   CORE FUNCTIONS — Orientation (rotation)
# #############################################################################

detect_animal_pole_robust <- function(spots) {
  early <- spots %>% filter(FRAME == min(FRAME, na.rm = TRUE))
  if (nrow(early) == 0) early <- spots
  c(mean(early$POSITION_X, na.rm = TRUE),
    mean(early$POSITION_Y, na.rm = TRUE),
    mean(early$POSITION_Z, na.rm = TRUE))
}

get_rotation_matrix <- function(u, v) {
  u <- u / sqrt(sum(u^2)); v <- v / sqrt(sum(v^2))
  if (sum(abs(u - v)) < 1e-6) return(diag(3))
  if (sum(abs(u + v)) < 1e-6) return(diag(c(1, -1, -1)))
  ax <- pracma::cross(u, v)
  ax <- ax / sqrt(sum(ax^2))
  ang <- acos(max(-1, min(1, sum(u * v))))
  K <- matrix(c(0, -ax[3], ax[2], ax[3], 0, -ax[1],
                -ax[2], ax[1], 0), 3, 3, byrow = TRUE)
  diag(3) + sin(ang) * K + (1 - cos(ang)) * (K %*% K)
}

apply_transform <- function(spots, ap, dorsal) {
  mid <- c(mean(spots$POSITION_X, na.rm = TRUE),
           mean(spots$POSITION_Y, na.rm = TRUE),
           mean(spots$POSITION_Z, na.rm = TRUE))
  mat <- as.matrix(spots[, c("POSITION_X", "POSITION_Y", "POSITION_Z")])
  mat_c <- sweep(mat, 2, mid)

  # Animal pole → +Y
  R1 <- get_rotation_matrix(ap - mid, c(0, 1, 0))
  mat_r1 <- t(R1 %*% t(mat_c))

  # Dorsal → +X (preserving +Y)
  d_r1 <- as.numeric(R1 %*% (dorsal - mid))
  d_xz <- c(d_r1[1], 0, d_r1[3])
  R2 <- if (sum(abs(d_xz)) < 1e-8) diag(3) else get_rotation_matrix(d_xz, c(1, 0, 0))
  mat_f <- t(R2 %*% t(mat_r1))

  spots$POSITION_X <- mat_f[, 1]
  spots$POSITION_Y <- mat_f[, 2]
  spots$POSITION_Z <- mat_f[, 3]
  list(data = spots, R1 = R1, R2 = R2, center = mid)
}

transform_tracks <- function(tracks_df, center, R1, R2) {
  if (is.null(tracks_df) || nrow(tracks_df) == 0) return(tracks_df)
  cols <- colnames(tracks_df)
  cand <- grep("_(?:POSITION_)?[XYZ]$", cols, perl = TRUE, value = TRUE)
  if (length(cand) == 0) return(tracks_df)
  stems <- unique(sub("_(?:POSITION_)?[XYZ]$", "", cand, perl = TRUE))

  for (stem in stems) {
    xn <- c(paste0(stem, "_POSITION_X"), paste0(stem, "_X"))
    yn <- c(paste0(stem, "_POSITION_Y"), paste0(stem, "_Y"))
    zn <- c(paste0(stem, "_POSITION_Z"), paste0(stem, "_Z"))
    xname <- xn[which(xn %in% cols)[1]]
    yname <- yn[which(yn %in% cols)[1]]
    zname <- zn[which(zn %in% cols)[1]]
    if (is.na(xname) || is.na(yname) || is.na(zname)) next

    mat <- tryCatch(
      apply(tracks_df[, c(xname, yname, zname)], 2,
            function(col) suppressWarnings(as.numeric(col))),
      error = function(e) NULL)
    if (is.null(mat) || !is.matrix(mat) || ncol(mat) != 3) next

    mat_f <- t(R2 %*% t(R1 %*% t(sweep(mat, 2, center))))
    tracks_df[[xname]] <- mat_f[, 1]
    tracks_df[[yname]] <- mat_f[, 2]
    tracks_df[[zname]] <- mat_f[, 3]
  }
  tracks_df
}

# #############################################################################
#            CORE FUNCTIONS — Sphere Fitting (cap-aware, robust)
# #############################################################################
#
# SPIM lightsheet data captures a THIN SPHERICAL CAP of the embryo:
#   Z extent ≈ 20–35% of X/Y extent.
# The sphere center is far outside the data volume.
# Standard algebraic sphere fits fail because the thin slab gives
# an ill-conditioned linear system.
#
# Strategy:
#   1. Identify the thin axis (Z for all typical SPIM data).
#   2. Grid the two long axes and extract the outer surface
#      (95th percentile of the thin-axis coordinate per grid cell).
#   3. Use the spherical-cap formula  R = (a² + h²) / (2h)  for the
#      initial guess (a = base radius, h = cap height).
#   4. Refine with Nelder-Mead + L-BFGS-B nonlinear least squares.
# #############################################################################

#' Algebraic sphere fit (linearised least squares).
#' Only used as a secondary fallback for roughly spherical data.
fit_sphere_algebraic <- function(x, y, z) {
  n <- length(x)
  A <- cbind(2 * x, 2 * y, 2 * z, rep(1, n))
  b <- x^2 + y^2 + z^2
  res <- tryCatch(qr.solve(A, b), error = function(e) rep(NA_real_, 4))
  if (any(is.na(res))) return(list(center = c(NA, NA, NA), radius = NA, rmse = Inf))
  cx <- res[1]; cy <- res[2]; cz <- res[3]
  R <- sqrt(max(0, res[4] + cx^2 + cy^2 + cz^2))
  dists <- sqrt((x - cx)^2 + (y - cy)^2 + (z - cz)^2)
  list(center = c(cx, cy, cz), radius = R,
       residuals = dists - R, rmse = sqrt(mean((dists - R)^2)))
}

#' Extract the outer embryo surface from a thin slab of nuclei.
#'
#' For each cell in an (N x N) grid of the two longest axes, keeps the
#' 95th-percentile of the thin axis (= outermost nuclei).
#' Returns an Mx3 matrix of surface points.
#'
#' @param x,y,z  numeric vectors (all nuclei)
#' @param grid_n number of bins per long axis
#' @param surface_quantile quantile to extract (0.95 = outer 5%)
extract_surface <- function(x, y, z, grid_n = 40, surface_quantile = 0.95) {
  extents <- c(diff(range(x)), diff(range(y)), diff(range(z)))
  thin <- which.min(extents)

  # Reorder so thin axis is column 3
  if (thin == 1) { u <- y; v <- z; w <- x }
  else if (thin == 2) { u <- x; v <- z; w <- y }
  else { u <- x; v <- y; w <- z }

  # Grid on (u, v)
  u_brk <- seq(min(u) - 1e-6, max(u) + 1e-6, length.out = grid_n + 1)
  v_brk <- seq(min(v) - 1e-6, max(v) + 1e-6, length.out = grid_n + 1)
  u_bin <- findInterval(u, u_brk, all.inside = TRUE)
  v_bin <- findInterval(v, v_brk, all.inside = TRUE)

  dt <- data.table::data.table(u = u, v = v, w = w, ub = u_bin, vb = v_bin)

  # Upper (max-w) surface: outer embryo cap
  surf_hi <- dt[, {
    q <- quantile(w, surface_quantile)
    sel <- w >= q
    .(u = mean(u[sel]), v = mean(v[sel]), w = mean(w[sel]))
  }, by = .(ub, vb)]

  # Lower (min-w) surface
  surf_lo <- dt[, {
    q <- quantile(w, 1 - surface_quantile)
    sel <- w <= q
    .(u = mean(u[sel]), v = mean(v[sel]), w = mean(w[sel]))
  }, by = .(ub, vb)]

  # Map back to (x, y, z) order
  to_xyz <- function(su, sv, sw) {
    if (thin == 1) cbind(sw, su, sv)
    else if (thin == 2) cbind(su, sw, sv)
    else cbind(su, sv, sw)
  }

  list(
    upper = to_xyz(surf_hi$u, surf_hi$v, surf_hi$w),
    lower = to_xyz(surf_lo$u, surf_lo$v, surf_lo$w),
    thin_axis = thin
  )
}

#' Fit a sphere using nonlinear least squares.
#' Minimises  sum( (||p_i - c|| - R)^2 ).
#'
#' @param pts  Nx3 matrix of surface points
#' @param cx0,cy0,cz0,R0  initial guesses
#' @return list(center, radius, rmse, convergence)
fit_sphere_nls <- function(pts, cx0, cy0, cz0, R0) {
  x <- pts[, 1]; y <- pts[, 2]; z <- pts[, 3]

  obj <- function(par) {
    d <- sqrt((x - par[1])^2 + (y - par[2])^2 + (z - par[3])^2)
    sum((d - par[4])^2)
  }

  p0 <- c(cx0, cy0, cz0, R0)

  # Nelder-Mead (gradient-free, robust)
  r1 <- tryCatch(
    optim(p0, obj, method = "Nelder-Mead",
          control = list(maxit = 50000, reltol = 1e-12)),
    error = function(e) list(value = Inf, par = p0, convergence = 99))

  # L-BFGS-B with R > 0 constraint
  r2 <- tryCatch(
    optim(p0, obj, method = "L-BFGS-B",
          lower = c(-Inf, -Inf, -Inf, 1),
          control = list(maxit = 50000)),
    error = function(e) list(value = Inf, par = p0, convergence = 99))

  best <- if (r2$value < r1$value) r2 else r1
  cx <- best$par[1]; cy <- best$par[2]; cz <- best$par[3]
  R <- abs(best$par[4])
  d <- sqrt((x - cx)^2 + (y - cy)^2 + (z - cz)^2)

  list(center = c(cx, cy, cz), radius = R,
       residuals = d - R,
       rmse = sqrt(mean((d - R)^2)),
       convergence = best$convergence)
}

#' Compute initial sphere guess from spherical cap geometry.
#'
#' For a cap of height h and base radius a:
#'   R = (a² + h²) / (2h)
#'   center_thin = top_of_cap - R   (for upper surface)
#'
#' @param pts  Nx3 surface points
#' @param thin_axis which axis (1=X, 2=Y, 3=Z) is thin
#' @param side "upper" or "lower"
cap_initial_guess <- function(pts, thin_axis, side = "upper") {
  x <- pts[, 1]; y <- pts[, 2]; z <- pts[, 3]

  # The two long axes centroid
  cx0 <- mean(x); cy0 <- mean(y); cz0 <- mean(z)

  # Lateral radius from centroid of long axes
  if (thin_axis == 1) {
    a <- max(sqrt((y - cy0)^2 + (z - cz0)^2))
    h <- diff(range(x)); h <- max(h, 1)
    R0 <- (a^2 + h^2) / (2 * h)
    if (side == "upper") cx0 <- max(x) - R0 else cx0 <- min(x) + R0
  } else if (thin_axis == 2) {
    a <- max(sqrt((x - cx0)^2 + (z - cz0)^2))
    h <- diff(range(y)); h <- max(h, 1)
    R0 <- (a^2 + h^2) / (2 * h)
    if (side == "upper") cy0 <- max(y) - R0 else cy0 <- min(y) + R0
  } else {
    a <- max(sqrt((x - cx0)^2 + (y - cy0)^2))
    h <- diff(range(z)); h <- max(h, 1)
    R0 <- (a^2 + h^2) / (2 * h)
    if (side == "upper") cz0 <- max(z) - R0 else cz0 <- min(z) + R0
  }

  list(cx = cx0, cy = cy0, cz = cz0, R = R0, a = a, h = h)
}

#' Main sphere fitting entry point — cap-aware, robust.
#'
#' Extracts the outer embryo surface from the thin-slab SPIM data,
#' then fits a sphere using NLS with a cap-geometry initial guess.
#' Tries both upper and lower surfaces; returns the better fit.
#'
#' @param points  Nx3 matrix (x, y, z) of ALL nuclei positions
#' @param shell_quantile  surface extraction quantile (0.5–0.99)
#' @param grid_n  grid resolution for surface extraction
#' @param max_points  subsample input for speed
fit_sphere_robust <- function(points, shell_quantile = 0.75, n_iter = 3,
                              max_points = 100000,
                              grid_n = 40) {
  x <- points[, 1]; y <- points[, 2]; z <- points[, 3]

  # Subsample for speed
  if (length(x) > max_points) {
    set.seed(42)
    idx <- sample(length(x), max_points)
    x <- x[idx]; y <- y[idx]; z <- z[idx]
  }

  # Map shell_quantile (originally "keep outer N%") to surface quantile
  # shell_quantile=0.75 → keep outer 25% → surface_quantile = 0.95
  # shell_quantile=0.5  → keep outer 50% → surface_quantile = 0.90
  surface_q <- 0.80 + shell_quantile * 0.15  # map [0.5, 0.95] → [0.875, 0.9425]
  surface_q <- min(0.99, max(0.80, surface_q))

  # 1. Extract surfaces
  surfs <- extract_surface(x, y, z, grid_n = grid_n,
                           surface_quantile = surface_q)

  # 2. Fit upper surface
  g_hi <- cap_initial_guess(surfs$upper, surfs$thin_axis, side = "upper")
  fit_hi <- fit_sphere_nls(surfs$upper, g_hi$cx, g_hi$cy, g_hi$cz, g_hi$R)

  # 3. Fit lower surface
  g_lo <- cap_initial_guess(surfs$lower, surfs$thin_axis, side = "lower")
  fit_lo <- fit_sphere_nls(surfs$lower, g_lo$cx, g_lo$cy, g_lo$cz, g_lo$R)

  # 4. Choose the better fit (lower RMSE)
  fit <- if (fit_hi$rmse <= fit_lo$rmse) fit_hi else fit_lo
  used_surface <- if (fit_hi$rmse <= fit_lo$rmse) "upper" else "lower"

  # 5. Compute coverage and stats on ALL points
  fit$coverage <- estimate_coverage(x, y, z, fit$center, fit$radius)
  d_all    <- sqrt((x - fit$center[1])^2 + (y - fit$center[2])^2 +
                   (z - fit$center[3])^2)
  near     <- abs(d_all - fit$radius) < fit$radius * 0.20
  fit$rmse_surface <- if (sum(near) > 10)
    sqrt(mean((d_all[near] - fit$radius)^2)) else fit$rmse
  fit$n_surface <- sum(near)
  fit$n_total   <- length(x)
  fit$surface_used  <- used_surface
  fit$cap_height    <- if (used_surface == "upper") g_hi$h else g_lo$h
  fit$cap_base_rad  <- if (used_surface == "upper") g_hi$a else g_lo$a

  fit
}

#' Estimate fraction of sphere surface covered by data.
#' Returns value in [0, 1]; <0.5 suggests partial embryo imaging.
estimate_coverage <- function(x, y, z, center, R) {
  if (is.na(R) || R <= 0) return(0)
  dx <- x - center[1]; dy <- y - center[2]; dz <- z - center[3]
  r <- sqrt(dx^2 + dy^2 + dz^2)
  on_surf <- r > R * 0.5 & r < R * 1.5
  if (sum(on_surf) < 10) return(0)

  theta <- acos(pmin(1, pmax(-1, dy[on_surf] / pmax(r[on_surf], 1e-10))))
  phi   <- atan2(dz[on_surf], dx[on_surf])
  tb <- as.integer(cut(theta, seq(0, pi, length.out = 13), include.lowest = TRUE))
  pb <- as.integer(cut(phi,   seq(-pi, pi, length.out = 25), include.lowest = TRUE))
  length(unique(paste(tb, pb))) / (12 * 24)
}

# #############################################################################
#      CORE FUNCTIONS — Spherical Coordinates, Margin, Wireframe
# #############################################################################

#' Add spherical coordinates to spots.
#' After orientation: +Y = animal pole, +X = dorsal.
#'   theta = 0° at animal pole, 180° at vegetal pole
#'   phi   = 0° at dorsal (+X)
#'   spherical_depth = R - r (0 at surface, positive inward)
compute_spherical_coords <- function(spots, center, R) {
  dx <- spots$POSITION_X - center[1]
  dy <- spots$POSITION_Y - center[2]
  dz <- spots$POSITION_Z - center[3]
  r  <- sqrt(dx^2 + dy^2 + dz^2)
  r_safe <- pmax(r, 1e-10)

  spots %>% mutate(
    RADIAL_DIST          = r,
    SPHERICAL_DEPTH      = R - r,
    SPHERICAL_DEPTH_NORM = (R - r) / R,
    THETA_DEG            = acos(pmin(1, pmax(-1, dy / r_safe))) * 180 / pi,
    PHI_DEG              = atan2(dz, dx) * 180 / pi
  )
}

#' Flag spots within the margin latitude band (and optionally depth range).
flag_margin <- function(spots, theta_min, theta_max,
                        depth_filter = FALSE, depth_min = 0, depth_max = Inf) {
  spots <- spots %>%
    mutate(IN_MARGIN = THETA_DEG >= theta_min & THETA_DEG <= theta_max)
  if (depth_filter) {
    spots <- spots %>%
      mutate(IN_MARGIN = IN_MARGIN &
               SPHERICAL_DEPTH >= depth_min & SPHERICAL_DEPTH <= depth_max)
  }
  spots
}

#' Sphere wireframe for plotly overlay (single trace with NA breaks).
make_wireframe <- function(center, radius, n_mer = 12, n_par = 8, n_pts = 60) {
  cx <- center[1]; cy <- center[2]; cz <- center[3]; R <- radius
  lines <- list()
  for (i in seq_len(n_mer)) {
    phi <- (i - 1) * 2 * pi / n_mer
    theta <- seq(0, pi, length.out = n_pts)
    lines[[length(lines) + 1]] <- tibble(
      x = cx + R * sin(theta) * cos(phi),
      y = cy + R * cos(theta),
      z = cz + R * sin(theta) * sin(phi))
  }
  for (j in seq_len(n_par)) {
    theta <- j * pi / (n_par + 1)
    phi <- seq(0, 2 * pi, length.out = n_pts)
    lines[[length(lines) + 1]] <- tibble(
      x = cx + R * sin(theta) * cos(phi),
      y = cy + R * cos(theta),
      z = cz + R * sin(theta) * sin(phi))
  }
  # Insert NA breaks between each line segment
  sep <- tibble(x = NA_real_, y = NA_real_, z = NA_real_)
  do.call(rbind, lapply(lines, function(df) rbind(df, sep)))
}

#' Margin band circles (top & bottom latitude) for plotly.
make_margin_band <- function(center, radius, theta_min_deg, theta_max_deg,
                             n_phi = 80) {
  cx <- center[1]; cy <- center[2]; cz <- center[3]; R <- radius
  phi <- seq(0, 2 * pi, length.out = n_phi)
  mk_circle <- function(th_deg) {
    th <- th_deg * pi / 180
    tibble(x = cx + R * sin(th) * cos(phi),
           y = cy + R * cos(th),
           z = cz + R * sin(th) * sin(phi))
  }
  sep <- tibble(x = NA_real_, y = NA_real_, z = NA_real_)
  rbind(mk_circle(theta_min_deg), sep, mk_circle(theta_max_deg))
}

#' Build track-line data for plotly (NA breaks between tracks).
make_track_lines <- function(spots_sub, max_tracks = 300) {
  tids <- unique(spots_sub$TRACK_ID)
  if (length(tids) > max_tracks) { set.seed(42); tids <- sample(tids, max_tracks) }
  d <- spots_sub %>%
    filter(TRACK_ID %in% tids) %>%
    arrange(TRACK_ID, FRAME) %>%
    select(TRACK_ID, POSITION_X, POSITION_Y, POSITION_Z, FRAME)
  sep <- tibble(TRACK_ID = NA_real_, POSITION_X = NA_real_,
                POSITION_Y = NA_real_, POSITION_Z = NA_real_, FRAME = NA_real_)
  d_list <- split(d, d$TRACK_ID)
  do.call(rbind, lapply(d_list, function(df) rbind(df, sep)))
}


# #############################################################################
#                                   UI
# #############################################################################

ui <- fluidPage(
  tags$head(tags$style(HTML("
    .sidebar-panel { overflow-y: auto; max-height: 90vh; }
    .section-header { margin-top: 8px; margin-bottom: 4px; }
  "))),

  titlePanel("Embryo Orientation, Sphere Fit & Margin Selection"),
  sidebarLayout(
    sidebarPanel(
      width = 3, class = "sidebar-panel",

      # ── 1. Upload ──────────────────────────────────────────
      h4("1. Upload Data", class = "section-header"),
      fileInput("spots_upload", "Spots CSV", accept = ".csv"),
      fileInput("tracks_upload", "Tracks CSV (optional)", accept = ".csv"),
      hr(),

      # ── 2. View mode ───────────────────────────────────────
      h4("2. View Mode", class = "section-header"),
      radioButtons("view_mode", NULL,
                   choices = c("All spots" = "all",
                               "Tracks" = "tracks",
                               "Nuclei in time" = "temporal"),
                   selected = "all", inline = FALSE),

      # Tracks mode options
      conditionalPanel(
        condition = "input.view_mode == 'tracks'",
        sliderInput("track_n_display", "Number of tracks",
                    min = 10, max = 2000, value = 300, step = 10),
        sliderInput("track_min_length", "Min track length (frames)",
                    min = 1, max = 100, value = 5, step = 1),
        selectInput("track_color", "Colour tracks by",
                    choices = c("Frame (start)" = "frame_start",
                                "Mean speed" = "mean_speed",
                                "Track duration" = "duration",
                                "Z position" = "z_mean"),
                    selected = "frame_start")
      ),

      # Nuclei-in-time mode options
      conditionalPanel(
        condition = "input.view_mode == 'temporal'",
        uiOutput("time_slider_ui"),
        checkboxInput("show_tails", "Show trailing tails", value = TRUE),
        sliderInput("track_tail", "Tail length (frames)",
                    min = 1, max = 50, value = 10, step = 1)
      ),
      hr(),

      # ── 3. Landmarks ───────────────────────────────────────
      h4("3. Pick Landmarks", class = "section-header"),
      p(style = "font-size:0.85em;",
        "Toggle a button, then click a dot on Tab 1."),
      uiOutput("ap_button"),
      textOutput("val_ap"),
      br(),
      uiOutput("dorsal_button"),
      textOutput("val_dorsal"),
      hr(),

      # ── Display ─────────────────────────────────────────────
      h4("Display", class = "section-header"),
      sliderInput("plot_n_spots", "Spots to display",
                  min = 500, max = 100000, value = 10000, step = 500),
      sliderInput("marker_size", "Marker size", min = 1, max = 10, value = 2),

      # Colour selector (context-dependent)
      conditionalPanel(
        condition = "input.main_tabs == 'oriented'",
        selectInput("oriented_color", "Colour by",
                    choices = c("Z (depth)" = "POSITION_Z",
                                "Dorsal–Ventral (X)" = "POSITION_X",
                                "Animal–Vegetal (Y)" = "POSITION_Y",
                                "Frame" = "FRAME",
                                "Spherical Depth" = "SPHERICAL_DEPTH",
                                "Latitude from AP" = "THETA_DEG"),
                    selected = "POSITION_Z")
      ),
      conditionalPanel(
        condition = "input.main_tabs == 'depth'",
        selectInput("depth_color", "Colour by",
                    choices = c("Spherical Depth" = "SPHERICAL_DEPTH",
                                "Normalised Depth" = "SPHERICAL_DEPTH_NORM",
                                "Latitude from AP" = "THETA_DEG",
                                "Longitude" = "PHI_DEG",
                                "Margin Status" = "margin_status"),
                    selected = "SPHERICAL_DEPTH")
      ),
      hr(),

      # ── 4. Orient ──────────────────────────────────────────
      h4("4. Orient Embryo", class = "section-header"),
      actionButton("btn_process", "Standardise Orientation",
                   class = "btn-success", width = "100%"),
      hr(),

      # ── 5. Sphere Fit ──────────────────────────────────────
      h4("5. Fit Sphere", class = "section-header"),
      sliderInput("shell_quantile", "Surface shell quantile",
                  min = 0.5, max = 0.95, value = 0.75, step = 0.05),
      helpText(style = "font-size:0.8em;",
               "Higher = stricter surface selection. Lower values tolerate ",
               "noisy / partial data."),
      actionButton("btn_fit_sphere", "Fit Sphere",
                   class = "btn-primary", width = "100%"),
      verbatimTextOutput("sphere_info"),
      hr(),

      # ── 6. Margin Selection ────────────────────────────────
      h4("6. Select Margin", class = "section-header"),
      helpText(style = "font-size:0.8em;",
               "Latitude: 0° = animal pole, 90° = equator, 180° = vegetal pole. ",
               "Ingression typically occurs in a band around the margin (60–130°)."),
      sliderInput("margin_theta", "Latitude band (degrees from AP)",
                  min = 0, max = 180, value = c(60, 120), step = 1),
      checkboxInput("margin_depth_filter", "Also filter by spherical depth",
                    value = FALSE),
      conditionalPanel(
        condition = "input.margin_depth_filter",
        sliderInput("margin_depth_range", "Spherical depth range (µm)",
                    min = -50, max = 200, value = c(0, 50), step = 1)
      ),
      verbatimTextOutput("margin_info"),
      hr(),

      # ── 7. Export ──────────────────────────────────────────
      h4("7. Export", class = "section-header"),
      downloadButton("download_data", "Download Enriched Data (ZIP)",
                     style = "width:100%;")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Explore Data", value = "raw",
                 plotlyOutput("plot_interactive", height = "700px")),
        tabPanel("Oriented & Sphere", value = "oriented",
                 plotlyOutput("plot_oriented", height = "700px")),
        tabPanel("Depth & Margin", value = "depth",
                 plotlyOutput("plot_depth_margin", height = "700px"))
      )
    )
  )
)


# #############################################################################
#                                  SERVER
# #############################################################################

server <- function(input, output, session) {

  # ═══════════════════ Reactive values ═══════════════════════════════════════
  vals <- reactiveValues(
    spots = NULL, tracks = NULL,
    animal_pole = NULL, dorsal_landmark = NULL,
    pick_mode = "none",
    oriented_spots = NULL, oriented_tracks = NULL,
    transform_params = NULL,
    sphere_fit = NULL,
    enriched_spots = NULL,   # oriented + spherical coords
    camera_raw = NULL, camera_ori = NULL, camera_dep = NULL
  )

  # ═══════════════════ 1. DATA LOADING ══════════════════════════════════════
  observe({
    req(input$spots_upload)
    tryCatch({
      tracks_path <- if (!is.null(input$tracks_upload)) input$tracks_upload$datapath else NULL
      data <- load_trackmate_data(input$spots_upload$datapath, tracks_path)
      vals$spots  <- data$spots
      vals$tracks <- data$tracks
      vals$animal_pole <- detect_animal_pole_robust(vals$spots)
      # Reset downstream
      vals$dorsal_landmark <- NULL
      vals$oriented_spots  <- NULL; vals$oriented_tracks <- NULL
      vals$sphere_fit <- NULL; vals$enriched_spots <- NULL

      n_s <- format(nrow(vals$spots), big.mark = ",")
      n_t <- if (!is.null(vals$tracks)) format(nrow(vals$tracks), big.mark = ",") else "none"
      n_f <- n_distinct(vals$spots$FRAME)
      showNotification(sprintf("Loaded %s spots, %s tracks, %d frames", n_s, n_t, n_f),
                       type = "message")
    }, error = function(e) {
      showNotification(paste("Load error:", e$message), type = "error")
    })
  })

  # ═══════════════════ 2. TIME SLIDER ═══════════════════════════════════════
  output$time_slider_ui <- renderUI({
    req(vals$spots)
    fmin <- as.integer(min(vals$spots$FRAME, na.rm = TRUE))
    fmax <- as.integer(max(vals$spots$FRAME, na.rm = TRUE))
    sliderInput("time_frame", "Frame",
                min = fmin, max = fmax, value = fmin, step = 1,
                animate = animationOptions(interval = 300, loop = TRUE))
  })

  # Debounced frame value for smooth temporal scrubbing
  time_frame_d <- debounce(reactive(input$time_frame), 200)

  # ═══════════════════ 3. LANDMARK PICKING ═════════════════════════════════
  output$ap_button <- renderUI({
    cls <- if (vals$pick_mode == "ap") "btn-warning" else "btn-primary"
    lbl <- if (vals$pick_mode == "ap") "Animal Pole (SELECTING)" else "Set Animal Pole"
    actionButton("btn_pick_ap", lbl, width = "100%", class = cls)
  })
  output$dorsal_button <- renderUI({
    cls <- if (vals$pick_mode == "dorsal") "btn-warning" else "btn-info"
    lbl <- if (vals$pick_mode == "dorsal") "Dorsal (SELECTING)" else "Set Dorsal"
    actionButton("btn_pick_dorsal", lbl, width = "100%", class = cls)
  })

  observeEvent(input$btn_pick_ap, {
    vals$pick_mode <- if (vals$pick_mode == "ap") "none" else "ap"
  })
  observeEvent(input$btn_pick_dorsal, {
    vals$pick_mode <- if (vals$pick_mode == "dorsal") "none" else "dorsal"
  })

  # Click handler — works on Tab 1 (source = "raw_plot")
  observeEvent(event_data("plotly_click", source = "raw_plot"), {
    cd <- event_data("plotly_click", source = "raw_plot")
    if (is.null(cd)) return()
    coords <- c(cd$x, cd$y, cd$z)
    if (vals$pick_mode == "ap") {
      vals$animal_pole <- coords; vals$pick_mode <- "none"
      showNotification("Animal Pole set")
    } else if (vals$pick_mode == "dorsal") {
      vals$dorsal_landmark <- coords; vals$pick_mode <- "none"
      showNotification("Dorsal landmark set")
    }
  })

  output$val_ap <- renderText({
    if (is.null(vals$animal_pole)) "Not set"
    else sprintf("AP: (%.1f, %.1f, %.1f)", vals$animal_pole[1],
                 vals$animal_pole[2], vals$animal_pole[3])
  })
  output$val_dorsal <- renderText({
    if (is.null(vals$dorsal_landmark)) "Not set"
    else sprintf("Dorsal: (%.1f, %.1f, %.1f)", vals$dorsal_landmark[1],
                 vals$dorsal_landmark[2], vals$dorsal_landmark[3])
  })

  # ─── Camera capture (raw plot) ───
  observeEvent(event_data("plotly_relayout", source = "raw_plot"), {
    rel <- event_data("plotly_relayout", source = "raw_plot")
    if (is.null(rel)) return()
    nm <- names(rel)
    ex <- rel[[ nm[grepl("camera.eye.x", nm)][1] ]]
    ey <- rel[[ nm[grepl("camera.eye.y", nm)][1] ]]
    ez <- rel[[ nm[grepl("camera.eye.z", nm)][1] ]]
    if (!is.null(ex) && !is.null(ey) && !is.null(ez))
      vals$camera_raw <- list(eye = list(x = ex, y = ey, z = ez))
  })

  # ═══════════════════ 4. ORIENTATION ══════════════════════════════════════
  observeEvent(input$btn_process, {
    req(vals$spots, vals$animal_pole, vals$dorsal_landmark)
    withProgress(message = "Orienting embryo...", {
      res <- apply_transform(vals$spots, vals$animal_pole, vals$dorsal_landmark)
      vals$oriented_spots  <- res$data
      vals$oriented_tracks <- transform_tracks(vals$tracks, res$center, res$R1, res$R2)
      vals$transform_params <- list(center = res$center, R1 = res$R1, R2 = res$R2)
      vals$sphere_fit <- NULL; vals$enriched_spots <- NULL
    })
    updateTabsetPanel(session, "main_tabs", selected = "oriented")
    showNotification("Orientation complete — now fit sphere.", type = "message")
  })

  # ═══════════════════ 5. SPHERE FITTING ═══════════════════════════════════
  observeEvent(input$btn_fit_sphere, {
    req(vals$oriented_spots)
    pts <- as.matrix(vals$oriented_spots[, c("POSITION_X", "POSITION_Y", "POSITION_Z")])
    fit <- withProgress(message = "Fitting sphere...", {
      fit_sphere_robust(pts, shell_quantile = input$shell_quantile)
    })

    vals$sphere_fit <- fit
    vals$enriched_spots <- compute_spherical_coords(
      vals$oriented_spots, fit$center, fit$radius)

    # Update depth slider range dynamically
    d_range <- range(vals$enriched_spots$SPHERICAL_DEPTH, na.rm = TRUE)
    updateSliderInput(session, "margin_depth_range",
                      min = floor(d_range[1]), max = ceiling(d_range[2]),
                      value = c(max(0, floor(d_range[1])),
                                ceiling(d_range[2] * 0.3)))

    if (!is.na(fit$coverage) && fit$coverage < 0.4) {
      showNotification(
        sprintf("Low sphere coverage (%.0f%%). Partial embryo detected — verify fit visually.",
                fit$coverage * 100),
        type = "warning", duration = 10)
    }
    showNotification(
      sprintf("Sphere: R=%.1f µm | RMSE=%.2f | Coverage=%.0f%%",
              fit$radius, fit$rmse, fit$coverage * 100),
      type = "message")
  })

  output$sphere_info <- renderText({
    sf <- vals$sphere_fit
    if (is.null(sf)) return("Orient embryo first, then fit sphere.")
    cap_info <- if (!is.null(sf$cap_height))
      sprintf("\nCap: h=%.0f µm, base_a=%.0f µm (%s surface)",
              sf$cap_height, sf$cap_base_rad,
              if (!is.null(sf$surface_used)) sf$surface_used else "?")
    else ""
    sprintf("R = %.1f µm | RMSE = %.2f µm\nCoverage = %.0f%% (%s)\nCentre: (%.1f, %.1f, %.1f)%s",
            sf$radius, sf$rmse, sf$coverage * 100,
            if (sf$coverage < 0.4) "PARTIAL — check visually"
            else if (sf$coverage < 0.7) "partial"
            else "good",
            sf$center[1], sf$center[2], sf$center[3],
            cap_info)
  })

  # ═══════════════════ 6. MARGIN FLAGGING ═════════════════════════════════
  # Reactive: enriched spots with current margin selection applied
  enriched_with_margin <- reactive({
    req(vals$enriched_spots)
    theta_rng <- input$margin_theta
    flag_margin(vals$enriched_spots, theta_rng[1], theta_rng[2],
                depth_filter  = isTRUE(input$margin_depth_filter),
                depth_min = if (!is.null(input$margin_depth_range))
                  input$margin_depth_range[1] else 0,
                depth_max = if (!is.null(input$margin_depth_range))
                  input$margin_depth_range[2] else Inf)
  })

  output$margin_info <- renderText({
    if (is.null(vals$enriched_spots)) return("Fit sphere first.")
    dat <- enriched_with_margin()
    n_m <- sum(dat$IN_MARGIN, na.rm = TRUE)
    n_t <- nrow(dat)
    theta_rng <- input$margin_theta
    sprintf("Margin: %s / %s spots (%.1f%%)\nLatitude band: %.0f°–%.0f° from AP",
            format(n_m, big.mark = ","), format(n_t, big.mark = ","),
            n_m / n_t * 100, theta_rng[1], theta_rng[2])
  })

  # ═══════════════════════════════════════════════════════════════════════════
  #                         TAB 1: RAW DATA & TIME
  # ═══════════════════════════════════════════════════════════════════════════

  output$plot_interactive <- renderPlotly({
    req(vals$spots)
    all_spots <- vals$spots
    msize <- max(1, input$marker_size)
    vmode <- input$view_mode  # "all", "tracks", "temporal"

    p <- plot_ly(source = "raw_plot")
    ttl <- ""  # title built per mode

    # ═════════════════════════════════════════════════════════════════════
    #  MODE 1: ALL SPOTS — static overview coloured by frame
    # ═════════════════════════════════════════════════════════════════════
    if (vmode == "all") {
      nshow <- min(nrow(all_spots), input$plot_n_spots)
      pd <- if (nrow(all_spots) > nshow) all_spots[sample(nrow(all_spots), nshow), ]
            else all_spots

      p <- p %>%
        add_trace(data = pd, x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
                  type = "scatter3d", mode = "markers",
                  marker = list(size = msize, color = ~FRAME,
                                colorscale = "Viridis", showscale = TRUE,
                                colorbar = list(title = "Frame")),
                  hoverinfo = "text",
                  text = ~paste("ID:", ID, "| Track:", TRACK_ID, "| F:", FRAME),
                  name = "Spots")

      ttl <- sprintf("All Spots — %s total", format(nrow(all_spots), big.mark = ","))

    # ═════════════════════════════════════════════════════════════════════
    #  MODE 2: TRACKS — full trajectories as 3D lines
    # ═════════════════════════════════════════════════════════════════════
    } else if (vmode == "tracks") {
      n_tracks <- if (!is.null(input$track_n_display)) input$track_n_display else 300
      min_len  <- if (!is.null(input$track_min_length)) input$track_min_length else 5
      col_by   <- if (!is.null(input$track_color)) input$track_color else "frame_start"

      # Compute per-track stats (only needs a few columns)
      track_stats <- all_spots %>%
        group_by(TRACK_ID) %>%
        summarise(n_frames    = n(),
                  frame_start = min(FRAME, na.rm = TRUE),
                  z_mean      = mean(POSITION_Z, na.rm = TRUE),
                  .groups = "drop") %>%
        filter(n_frames >= min_len)

      # Optionally add mean speed if colouring by it
      if (col_by == "mean_speed" || col_by == "duration") {
        track_vel <- all_spots %>%
          arrange(TRACK_ID, FRAME) %>%
          group_by(TRACK_ID) %>%
          mutate(disp = sqrt((POSITION_X - lag(POSITION_X))^2 +
                             (POSITION_Y - lag(POSITION_Y))^2 +
                             (POSITION_Z - lag(POSITION_Z))^2)) %>%
          summarise(mean_speed = mean(disp, na.rm = TRUE),
                    duration   = max(FRAME) - min(FRAME),
                    .groups = "drop")
        track_stats <- track_stats %>% left_join(track_vel, by = "TRACK_ID")
      }

      # Sample tracks
      if (nrow(track_stats) > n_tracks) {
        set.seed(42)
        track_stats <- track_stats[sample(nrow(track_stats), n_tracks), ]
      }
      tids <- track_stats$TRACK_ID

      # Get ordered spots for selected tracks
      track_spots <- all_spots %>%
        filter(TRACK_ID %in% tids) %>%
        left_join(track_stats %>% select(TRACK_ID, any_of(col_by)),
                  by = "TRACK_ID") %>%
        arrange(TRACK_ID, FRAME)

      color_val <- track_spots[[col_by]]
      cname <- switch(col_by,
                      frame_start = "Start frame", mean_speed = "Speed (µm/step)",
                      duration = "Duration (frames)", z_mean = "Mean Z (µm)", col_by)

      # Build line segments with NA breaks between tracks
      tl <- make_track_lines(track_spots, max_tracks = n_tracks + 1)
      # Propagate colour from spots to line data
      tl_col <- tl %>%
        left_join(track_stats %>% select(TRACK_ID, any_of(col_by)), by = "TRACK_ID")

      p <- p %>%
        add_trace(data = tl_col, x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
                  type = "scatter3d", mode = "lines",
                  line = list(color = tl_col[[col_by]],
                              colorscale = "Plasma", width = 2,
                              showscale = TRUE,
                              colorbar = list(title = cname)),
                  hoverinfo = "none", showlegend = FALSE)

      # Add small start-point markers so user can click for landmark picking
      start_pts <- track_spots %>%
        group_by(TRACK_ID) %>% slice_min(FRAME, n = 1) %>% ungroup()
      nshow_pts <- min(nrow(start_pts), input$plot_n_spots)
      if (nrow(start_pts) > nshow_pts)
        start_pts <- start_pts[sample(nrow(start_pts), nshow_pts), ]

      p <- p %>%
        add_trace(data = start_pts,
                  x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
                  type = "scatter3d", mode = "markers",
                  marker = list(size = max(1, msize - 1), color = "grey60", opacity = 0.4),
                  hoverinfo = "text",
                  text = ~paste("Track:", TRACK_ID, "| F:", FRAME),
                  name = "Track starts")

      ttl <- sprintf("%d tracks (≥%d frames) — colour: %s",
                     length(tids), min_len, cname)

    # ═════════════════════════════════════════════════════════════════════
    #  MODE 3: NUCLEI IN TIME — scrub through frames, watch movement
    # ═════════════════════════════════════════════════════════════════════
    } else { # vmode == "temporal"
      fv <- time_frame_d()
      req(fv)
      current <- all_spots %>% filter(FRAME == fv)

      # Subsample current frame
      nshow <- min(nrow(current), input$plot_n_spots)
      pd <- if (nrow(current) > nshow) current[sample(nrow(current), nshow), ]
            else current

      p <- p %>%
        add_trace(data = pd, x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
                  type = "scatter3d", mode = "markers",
                  marker = list(size = msize, color = "#2166AC"),
                  hoverinfo = "text",
                  text = ~paste("ID:", ID, "| Track:", TRACK_ID),
                  name = sprintf("Frame %d", fv))

      # Trailing tails
      if (isTRUE(input$show_tails)) {
        tl_len <- if (!is.null(input$track_tail)) input$track_tail else 5
        f_min <- max(min(all_spots$FRAME, na.rm = TRUE), fv - tl_len)
        tids <- unique(current$TRACK_ID)
        tail_data <- all_spots %>%
          filter(TRACK_ID %in% tids & FRAME >= f_min & FRAME <= fv)
        if (nrow(tail_data) > 0) {
          tl <- make_track_lines(tail_data, max_tracks = 500)
          p <- p %>%
            add_trace(data = tl,
                      x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
                      type = "scatter3d", mode = "lines",
                      line = list(color = "rgba(255,200,0,0.5)", width = 2),
                      hoverinfo = "none", name = "Tails")
        }
      }

      ttl <- sprintf("Frame %d / %d (%.1f min) — %d nuclei",
                     fv, max(all_spots$FRAME, na.rm = TRUE),
                     fv * 0.5, nrow(current))
    }

    # ── Landmarks (all modes) ───────────────────────────────────────────
    if (!is.null(vals$animal_pole)) {
      p <- p %>% add_trace(x = vals$animal_pole[1], y = vals$animal_pole[2],
                           z = vals$animal_pole[3],
                           type = "scatter3d", mode = "markers",
                           marker = list(size = msize + 6, color = "red",
                                         symbol = "diamond"),
                           name = "Animal Pole")
    }
    if (!is.null(vals$dorsal_landmark)) {
      p <- p %>% add_trace(x = vals$dorsal_landmark[1], y = vals$dorsal_landmark[2],
                           z = vals$dorsal_landmark[3],
                           type = "scatter3d", mode = "markers",
                           marker = list(size = msize + 6, color = "cyan",
                                         symbol = "square"),
                           name = "Dorsal")
    }

    # ── Title override for picking mode ─────────────────────────────────
    pick <- isolate(vals$pick_mode)
    if (pick == "ap") ttl <- "CLICK A SPOT TO SET ANIMAL POLE"
    else if (pick == "dorsal") ttl <- "CLICK A SPOT TO SET DORSAL"

    scene <- list(aspectmode = "data")
    if (!is.null(vals$camera_raw)) scene$camera <- vals$camera_raw

    p %>% layout(title = list(text = ttl), scene = scene) %>%
      event_register("plotly_click")
  })

  # ═══════════════════════════════════════════════════════════════════════════
  #                     TAB 2: ORIENTED & SPHERE
  # ═══════════════════════════════════════════════════════════════════════════

  output$plot_oriented <- renderPlotly({
    req(vals$oriented_spots)
    dat <- if (!is.null(vals$enriched_spots)) vals$enriched_spots
           else vals$oriented_spots
    msize <- max(1, input$marker_size)
    nshow <- min(nrow(dat), input$plot_n_spots)
    if (nrow(dat) > nshow) dat <- dat[sample(nrow(dat), nshow), ]

    # Colour
    cchoice <- if (!is.null(input$oriented_color)) input$oriented_color else "POSITION_Z"
    if (!cchoice %in% names(dat)) cchoice <- "POSITION_Z"
    cval <- dat[[cchoice]]
    cname <- switch(cchoice,
                    POSITION_Z = "Z", POSITION_X = "X (D-V)", POSITION_Y = "Y (AP-VP)",
                    FRAME = "Frame", SPHERICAL_DEPTH = "Sph. Depth (µm)",
                    THETA_DEG = "Latitude (°)", cchoice)

    p <- plot_ly(source = "oriented_plot") %>%
      add_trace(data = dat, x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
                type = "scatter3d", mode = "markers",
                marker = list(size = msize, color = cval,
                              colorscale = "Plasma", showscale = TRUE,
                              colorbar = list(title = cname)),
                hoverinfo = "text",
                text = ~paste("Track:", TRACK_ID, "| F:", FRAME),
                name = "Spots")

    # Sphere wireframe
    sf <- vals$sphere_fit
    if (!is.null(sf)) {
      wf <- make_wireframe(sf$center, sf$radius)
      p <- p %>%
        add_trace(data = wf, x = ~x, y = ~y, z = ~z,
                  type = "scatter3d", mode = "lines",
                  line = list(color = "rgba(120,120,120,0.25)", width = 1),
                  hoverinfo = "none", name = "Sphere fit") %>%
        add_trace(x = sf$center[1], y = sf$center[2], z = sf$center[3],
                  type = "scatter3d", mode = "markers",
                  marker = list(size = 6, color = "orange", symbol = "cross"),
                  name = "Sphere centre")

      # Margin band
      theta_rng <- input$margin_theta
      mb <- make_margin_band(sf$center, sf$radius, theta_rng[1], theta_rng[2])
      p <- p %>%
        add_trace(data = mb, x = ~x, y = ~y, z = ~z,
                  type = "scatter3d", mode = "lines",
                  line = list(color = "#E41A1C", width = 3),
                  hoverinfo = "none", name = "Margin band")
    }

    sub <- if (!is.null(sf))
      sprintf("R=%.0f µm | RMSE=%.1f | Coverage=%.0f%%",
              sf$radius, sf$rmse, sf$coverage * 100)
    else "Orient done. Fit sphere (step 5) to add wireframe."

    scene <- list(aspectmode = "data",
                  xaxis = list(title = "Dorsal (+X)"),
                  yaxis = list(title = "AP (+Y)"),
                  zaxis = list(title = "Z"))

    p %>% layout(title = list(text = paste0("Oriented<br><sup>", sub, "</sup>")),
                 scene = scene)
  })

  # ═══════════════════════════════════════════════════════════════════════════
  #                      TAB 3: DEPTH & MARGIN
  # ═══════════════════════════════════════════════════════════════════════════

  output$plot_depth_margin <- renderPlotly({
    req(vals$enriched_spots)
    sf <- vals$sphere_fit
    theta_rng <- input$margin_theta
    msize <- max(1, input$marker_size)

    # Apply margin flag
    full <- enriched_with_margin()
    n_tot <- nrow(full)
    n_margin <- sum(full$IN_MARGIN, na.rm = TRUE)

    # --- Subsample (oversample margin for visibility) ---
    nshow <- min(n_tot, input$plot_n_spots)
    if (n_tot > nshow) {
      if (n_margin > 0 && n_margin < nshow * 0.5) {
        m_show <- min(n_margin, round(nshow * 0.5))
        o_show <- nshow - m_show
        m_idx <- sample(which(full$IN_MARGIN), m_show)
        o_idx <- sample(which(!full$IN_MARGIN), min(o_show, sum(!full$IN_MARGIN)))
        dat <- full[c(m_idx, o_idx), ]
      } else {
        dat <- full[sample(n_tot, nshow), ]
      }
    } else {
      dat <- full
    }

    # Colour mode
    cmode <- if (!is.null(input$depth_color)) input$depth_color else "SPHERICAL_DEPTH"

    if (cmode == "margin_status") {
      # ── Binary colouring: margin vs non-margin ──
      p <- plot_ly(source = "depth_plot") %>%
        add_trace(data = dat %>% filter(!IN_MARGIN),
                  x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
                  type = "scatter3d", mode = "markers",
                  marker = list(size = msize, color = "grey75", opacity = 0.3),
                  name = "Non-margin", hoverinfo = "text",
                  text = ~sprintf("θ=%.0f° d=%.1fµm", THETA_DEG, SPHERICAL_DEPTH)) %>%
        add_trace(data = dat %>% filter(IN_MARGIN),
                  x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
                  type = "scatter3d", mode = "markers",
                  marker = list(size = msize + 1, color = "#E41A1C"),
                  name = "Margin", hoverinfo = "text",
                  text = ~sprintf("θ=%.0f° d=%.1fµm", THETA_DEG, SPHERICAL_DEPTH))
    } else {
      # ── Continuous colouring ──
      cval <- dat[[cmode]]
      cname <- switch(cmode,
                      SPHERICAL_DEPTH = "Depth (µm)", SPHERICAL_DEPTH_NORM = "Norm. Depth",
                      THETA_DEG = "Latitude (°)", PHI_DEG = "Longitude (°)",
                      RADIAL_DIST = "Radial Dist (µm)", cmode)
      p <- plot_ly(source = "depth_plot") %>%
        add_trace(data = dat, x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
                  type = "scatter3d", mode = "markers",
                  marker = list(size = msize, color = cval,
                                colorscale = "Plasma", showscale = TRUE,
                                colorbar = list(title = cname)),
                  hoverinfo = "text",
                  text = ~sprintf("θ=%.0f° d=%.1fµm F=%d",
                                  THETA_DEG, SPHERICAL_DEPTH, FRAME),
                  name = "Spots")

      # Highlight margin boundary with rings around margin spots
      if (any(dat$IN_MARGIN)) {
        mp <- dat %>% filter(IN_MARGIN)
        p <- p %>%
          add_trace(data = mp, x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
                    type = "scatter3d", mode = "markers",
                    marker = list(size = msize + 2,
                                  color = "rgba(0,0,0,0)",
                                  line = list(color = "#E41A1C", width = 2)),
                    name = "Margin", hoverinfo = "none")
      }
    }

    # ── Sphere wireframe ──
    if (!is.null(sf)) {
      wf <- make_wireframe(sf$center, sf$radius)
      p <- p %>%
        add_trace(data = wf, x = ~x, y = ~y, z = ~z,
                  type = "scatter3d", mode = "lines",
                  line = list(color = "rgba(100,100,100,0.2)", width = 1),
                  hoverinfo = "none", name = "Sphere fit", showlegend = TRUE)

      mb <- make_margin_band(sf$center, sf$radius, theta_rng[1], theta_rng[2])
      p <- p %>%
        add_trace(data = mb, x = ~x, y = ~y, z = ~z,
                  type = "scatter3d", mode = "lines",
                  line = list(color = "#E41A1C", width = 3),
                  hoverinfo = "none", name = "Margin band")
    }

    scene <- list(aspectmode = "data",
                  xaxis = list(title = "Dorsal (+X)"),
                  yaxis = list(title = "AP (+Y)"),
                  zaxis = list(title = "Z"))

    sub <- sprintf("%s margin / %s total (%.1f%%) | θ ∈ [%.0f°, %.0f°]",
                   format(n_margin, big.mark = ","),
                   format(n_tot, big.mark = ","),
                   n_margin / n_tot * 100,
                   theta_rng[1], theta_rng[2])

    p %>% layout(
      title = list(text = paste0("Spherical Depth & Margin<br><sup>", sub, "</sup>")),
      scene = scene)
  })

  # ═══════════════════ 7. EXPORT ═══════════════════════════════════════════
  output$download_data <- downloadHandler(
    filename = function() paste0("oriented_data_", Sys.Date(), ".zip"),
    content = function(file) {
      tmpdir <- tempdir()
      files_to_zip <- character(0)

      # Best available spots
      spots_out <- if (!is.null(vals$enriched_spots)) {
        enriched_with_margin()
      } else if (!is.null(vals$oriented_spots)) {
        vals$oriented_spots
      } else {
        NULL
      }

      if (!is.null(spots_out)) {
        write_csv(spots_out, file.path(tmpdir, "oriented_spots.csv"))
        files_to_zip <- c(files_to_zip, "oriented_spots.csv")
      }
      if (!is.null(vals$oriented_tracks)) {
        write_csv(vals$oriented_tracks, file.path(tmpdir, "oriented_tracks.csv"))
        files_to_zip <- c(files_to_zip, "oriented_tracks.csv")
      }

      # Sphere parameters
      sf <- vals$sphere_fit
      if (!is.null(sf)) {
        theta_rng <- input$margin_theta
        params <- tibble(
          parameter = c("center_x", "center_y", "center_z",
                        "radius", "rmse", "coverage",
                        "margin_theta_min", "margin_theta_max",
                        "margin_depth_filter",
                        "margin_depth_min", "margin_depth_max"),
          value = c(sf$center[1], sf$center[2], sf$center[3],
                    sf$radius, sf$rmse, sf$coverage,
                    theta_rng[1], theta_rng[2],
                    as.numeric(isTRUE(input$margin_depth_filter)),
                    if (!is.null(input$margin_depth_range)) input$margin_depth_range[1] else NA,
                    if (!is.null(input$margin_depth_range)) input$margin_depth_range[2] else NA)
        )
        write_csv(params, file.path(tmpdir, "sphere_params.csv"))
        files_to_zip <- c(files_to_zip, "sphere_params.csv")
      }

      req(length(files_to_zip) > 0)
      oldwd <- setwd(tmpdir); on.exit(setwd(oldwd))
      utils::zip(zipfile = file, files = files_to_zip)
    }
  )
}

# #############################################################################
shinyApp(ui, server)
