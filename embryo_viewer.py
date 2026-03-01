# =============================================================================
# napari-based 4D Embryo Viewer — Orientation, Sphere Fit & Ingression Analysis
# =============================================================================
#
# PURPOSE:
#   Interactive 4D (3D + time) visualisation and annotation of light-sheet
#   microscopy nuclear-tracking data (TrackMate exports) for medaka and
#   zebrafish gastrulation analysis.
#
# REPLACES:  R/Shiny + plotly orientation app (too slow for frame-by-frame
#            scrubbing of millions of spots).
#
# FEATURES:
#   1. GPU-accelerated 4D point cloud (millions of nuclei, smooth scrubbing)
#   2. Progressive track formation as time advances
#   3. Interactive landmark picking (Animal Pole, Dorsal)
#   4. Robust sphere fitting for thin SPIM caps
#   5. Spherical coordinate system (depth, latitude, longitude)
#   6. ROI rectangle selection (draw in 2D, selects through full 3D volume)
#   7. Ingression analysis: radial velocity, ingression flagging & colouring
#   8. Export enriched CSVs → feeds R analysis pipeline
#
# USAGE:
#   python embryo_viewer.py medaka_25082025_combined_spots.csv  \
#                           --tracks medaka_25082025_combined_tracks.csv
#
# INSTALLATION:
#   pip install "napari[all]" scipy pandas
#
# OUTPUT (via Export button):
#   oriented_spots.csv  — all spots with spherical coords + IN_ROI flag
#   sphere_params.csv   — fitted sphere centre, radius
#
# PIPELINE:
#   → Step 1: THIS SCRIPT
#     Step 2: trackmate_filter_and_validate.R  (R)
#     Step 3: trackmate_analysis.R             (R)
#
# =============================================================================

from __future__ import annotations

import argparse
import os
import sys
import threading
from pathlib import Path

# ── Wayland / DPI fixes for Fedora GNOME ──────────────────────────────────
# Must be set before any Qt import
if os.environ.get("XDG_SESSION_TYPE") == "wayland":
    os.environ.setdefault("QT_QPA_PLATFORM", "wayland")
os.environ.setdefault("QT_AUTO_SCREEN_SCALE_FACTOR", "1")

import numpy as np
import pandas as pd
from scipy.optimize import minimize

# ---------------------------------------------------------------------------
#  napari import (deferred so --help works without Qt)
# ---------------------------------------------------------------------------
_napari = None
_magicgui = None


def _import_napari():
    global _napari, _magicgui
    if _napari is None:
        import napari as _nap
        import magicgui as _mg

        _napari = _nap
        _magicgui = _mg
    return _napari, _magicgui


# ############################################################################
#                            DATA LOADING
# ############################################################################


def load_trackmate_csv(path: str | Path) -> pd.DataFrame:
    """Load a TrackMate spots/tracks CSV, auto-skipping metadata rows."""
    path = Path(path)
    df = pd.read_csv(path, low_memory=False)

    # TrackMate CSVs have 3 metadata rows after the header.
    # Detect by checking if the first value of POSITION_X (or similar
    # numeric column) is non-numeric.
    num_cols = [
        c
        for c in df.columns
        if "POSITION" in c or c in ("FRAME", "QUALITY", "RADIUS", "ID", "TRACK_ID")
    ]
    if len(num_cols) > 0:
        test_col = num_cols[0]
        try:
            float(df[test_col].iloc[0])
        except (ValueError, TypeError):
            # Drop first 3 rows (descriptions + units)
            df = df.iloc[3:].reset_index(drop=True)

    # Convert numeric columns
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    return df


def load_spots(path: str | Path) -> pd.DataFrame:
    """Load spots CSV and return cleaned DataFrame."""
    df = load_trackmate_csv(path)
    required = ["POSITION_X", "POSITION_Y", "POSITION_Z", "FRAME"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    df = df.dropna(subset=required)
    # Ensure integer frame
    df["FRAME"] = df["FRAME"].astype(int)
    return df


# ############################################################################
#                         SPHERE FITTING (cap-aware)
# ############################################################################


def extract_surface(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    grid_n: int = 40,
    quantile: float = 0.95,
) -> dict:
    """Extract outer and inner surface points from a SPIM thin slab.

    Returns dict with 'upper' and 'lower' surface point arrays (Nx3)
    and the identity of the thin axis.
    """
    extents = np.array([np.ptp(x), np.ptp(y), np.ptp(z)])
    thin_axis = int(np.argmin(extents))

    # Re-order so thin axis is last
    axes = [x, y, z]
    order = [i for i in range(3) if i != thin_axis] + [thin_axis]
    u, v, w = axes[order[0]], axes[order[1]], axes[order[2]]

    # Grid on (u, v)
    u_edges = np.linspace(u.min() - 1e-6, u.max() + 1e-6, grid_n + 1)
    v_edges = np.linspace(v.min() - 1e-6, v.max() + 1e-6, grid_n + 1)
    u_bin = np.digitize(u, u_edges) - 1
    v_bin = np.digitize(v, v_edges) - 1
    u_bin = np.clip(u_bin, 0, grid_n - 1)
    v_bin = np.clip(v_bin, 0, grid_n - 1)
    cell_id = u_bin * grid_n + v_bin

    upper_pts = []
    lower_pts = []
    for cid in np.unique(cell_id):
        mask = cell_id == cid
        if mask.sum() < 3:
            continue
        um, vm, wm = u[mask], v[mask], w[mask]
        # Upper surface: high-w
        q_hi = np.quantile(wm, quantile)
        sel_hi = wm >= q_hi
        upper_pts.append([um[sel_hi].mean(), vm[sel_hi].mean(), wm[sel_hi].mean()])
        # Lower surface: low-w
        q_lo = np.quantile(wm, 1 - quantile)
        sel_lo = wm <= q_lo
        lower_pts.append([um[sel_lo].mean(), vm[sel_lo].mean(), wm[sel_lo].mean()])

    upper = np.array(upper_pts)
    lower = np.array(lower_pts)

    # Map back to (x, y, z) order
    inv = [0, 0, 0]
    for new_i, orig_i in enumerate(order):
        inv[orig_i] = new_i

    def reorder(arr):
        return arr[:, inv]

    return {
        "upper": reorder(upper),
        "lower": reorder(lower),
        "thin_axis": thin_axis,
        "order": order,
    }


def _cap_initial_guess(pts: np.ndarray, thin_axis: int, side: str = "upper") -> dict:
    """Spherical-cap geometry initial guess for NLS.

    R = (a² + h²) / (2h)  where  a = base radius, h = cap height.
    """
    cx, cy, cz = pts.mean(axis=0)
    long_axes = [i for i in range(3) if i != thin_axis]
    center_lat = np.array([cx, cy, cz])

    # Lateral radius
    a = np.sqrt(
        (pts[:, long_axes[0]] - center_lat[long_axes[0]]) ** 2
        + (pts[:, long_axes[1]] - center_lat[long_axes[1]]) ** 2
    ).max()

    # Cap height along thin axis
    h = np.ptp(pts[:, thin_axis])
    h = max(h, 1.0)

    R = (a**2 + h**2) / (2 * h)

    if side == "upper":
        center_lat[thin_axis] = pts[:, thin_axis].max() - R
    else:
        center_lat[thin_axis] = pts[:, thin_axis].min() + R

    return {"center": center_lat, "R": R, "a": a, "h": h}


def _fit_sphere_nls(pts: np.ndarray, cx0, cy0, cz0, R0) -> dict:
    """Fit sphere by minimising sum((||p - c|| - R)^2)."""

    def obj(par):
        d = np.sqrt(
            (pts[:, 0] - par[0]) ** 2
            + (pts[:, 1] - par[1]) ** 2
            + (pts[:, 2] - par[2]) ** 2
        )
        return np.sum((d - par[3]) ** 2)

    p0 = np.array([cx0, cy0, cz0, R0], dtype=float)

    r1 = minimize(
        obj,
        p0,
        method="Nelder-Mead",
        options={"maxiter": 80000, "xatol": 1e-6, "fatol": 1e-10},
    )
    r2 = minimize(
        obj,
        p0,
        method="L-BFGS-B",
        bounds=[(None, None)] * 3 + [(1.0, None)],
        options={"maxiter": 50000},
    )

    best = r2 if r2.fun < r1.fun else r1
    cx, cy, cz, R = best.x
    R = abs(R)
    d = np.sqrt((pts[:, 0] - cx) ** 2 + (pts[:, 1] - cy) ** 2 + (pts[:, 2] - cz) ** 2)
    return {
        "center": np.array([cx, cy, cz]),
        "radius": R,
        "rmse": float(np.sqrt(np.mean((d - R) ** 2))),
        "residuals": d - R,
    }


def fit_sphere(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    max_points: int = 100_000,
    grid_n: int = 40,
    surface_quantile: float = 0.9125,
) -> dict:
    """Robust sphere fitting for thin SPIM caps.

    Parameters
    ----------
    x, y, z : position arrays (all nuclei)
    max_points : subsample for speed
    grid_n : grid resolution for surface extraction
    surface_quantile : top-% of thin axis to keep per grid cell

    Returns
    -------
    dict with center (3,), radius, rmse, coverage, surface_used, cap_height
    """
    if len(x) > max_points:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(x), max_points, replace=False)
        x, y, z = x[idx], y[idx], z[idx]

    surfs = extract_surface(x, y, z, grid_n=grid_n, quantile=surface_quantile)

    # Fit upper surface
    g_hi = _cap_initial_guess(surfs["upper"], surfs["thin_axis"], "upper")
    fit_hi = _fit_sphere_nls(surfs["upper"], *g_hi["center"], g_hi["R"])

    # Fit lower surface
    g_lo = _cap_initial_guess(surfs["lower"], surfs["thin_axis"], "lower")
    fit_lo = _fit_sphere_nls(surfs["lower"], *g_lo["center"], g_lo["R"])

    # Pick the better fit
    if fit_hi["rmse"] <= fit_lo["rmse"]:
        fit, side, guess = fit_hi, "upper", g_hi
    else:
        fit, side, guess = fit_lo, "lower", g_lo

    # Coverage estimate (fraction of sphere-surface bins hit)
    d = np.sqrt(
        (x - fit["center"][0]) ** 2
        + (y - fit["center"][1]) ** 2
        + (z - fit["center"][2]) ** 2
    )
    R = fit["radius"]
    on_surf = (d > R * 0.5) & (d < R * 1.5)
    coverage = 0.0
    if on_surf.sum() > 10:
        dx = x[on_surf] - fit["center"][0]
        dy = y[on_surf] - fit["center"][1]
        dz = z[on_surf] - fit["center"][2]
        r = np.sqrt(dx**2 + dy**2 + dz**2)
        r = np.maximum(r, 1e-10)
        theta = np.arccos(np.clip(-dy / r, -1, 1))  # AP at -Y
        phi = np.arctan2(dz, dx)
        tb = np.digitize(theta, np.linspace(0, np.pi, 13))
        pb = np.digitize(phi, np.linspace(-np.pi, np.pi, 25))
        bins = set(zip(tb.tolist(), pb.tolist()))
        coverage = len(bins) / (12 * 24)

    fit["coverage"] = coverage
    fit["surface_used"] = side
    fit["cap_height"] = guess["h"]
    fit["cap_base_radius"] = guess["a"]
    return fit


# ############################################################################
#                     ORIENTATION & SPHERICAL COORDS
# ############################################################################


def rotation_matrix(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Rotation matrix that maps unit vector u to unit vector v."""
    u = u / np.linalg.norm(u)
    v = v / np.linalg.norm(v)
    if np.allclose(u, v):
        return np.eye(3)
    if np.allclose(u, -v):
        return np.diag([1.0, -1.0, -1.0])
    ax = np.cross(u, v)
    ax = ax / np.linalg.norm(ax)
    ang = np.arccos(np.clip(u @ v, -1, 1))
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


def orient_embryo(
    pos: np.ndarray, ap: np.ndarray, dorsal: np.ndarray
) -> tuple[np.ndarray, dict]:
    """Rotate embryo so AP → -Y (top in napari 3D view), Dorsal → +X.

    Parameters
    ----------
    pos : (N, 3) array of positions
    ap : (3,) animal pole coordinate
    dorsal : (3,) dorsal landmark

    Returns
    -------
    (N, 3) oriented positions, dict of transform params
    """
    mid = pos.mean(axis=0)
    pos_c = pos - mid

    # Step 1: AP → -Y (napari +Y points downward in 3D, so -Y = top)
    R1 = rotation_matrix(ap - mid, np.array([0.0, -1.0, 0.0]))
    pos_r1 = (R1 @ pos_c.T).T

    # Step 2: Dorsal → +X (preserve +Y)
    d_r1 = R1 @ (dorsal - mid)
    d_xz = np.array([d_r1[0], 0.0, d_r1[2]])
    if np.linalg.norm(d_xz) < 1e-8:
        R2 = np.eye(3)
    else:
        R2 = rotation_matrix(d_xz, np.array([1.0, 0.0, 0.0]))

    pos_f = (R2 @ pos_r1.T).T
    return pos_f, {"center": mid, "R1": R1, "R2": R2}


def compute_spherical_coords(pos: np.ndarray, center: np.ndarray, R: float) -> dict:
    """Spherical coords relative to fitted sphere.

    Returns dict with keys:
      RADIAL_DIST, SPHERICAL_DEPTH, SPHERICAL_DEPTH_NORM,
      THETA_DEG (0=AP, 90=equator, 180=VP),
      PHI_DEG (0=dorsal +X direction)
    """
    d = pos - center
    r = np.sqrt(np.sum(d**2, axis=1))
    r_safe = np.maximum(r, 1e-10)
    return {
        "RADIAL_DIST": r,
        "SPHERICAL_DEPTH": R - r,
        "SPHERICAL_DEPTH_NORM": (R - r) / R,
        "THETA_DEG": np.degrees(np.arccos(np.clip(-d[:, 1] / r_safe, -1, 1))),
        "PHI_DEG": np.degrees(np.arctan2(d[:, 2], d[:, 0])),
    }


# ############################################################################
#                         SPHERE MESH (for napari)
# ############################################################################


def make_sphere_cap_mesh(
    center, radius, data_points, padding_deg=10, n_lat=30, n_lon=60
):
    """Create mesh for the sphere cap overlapping with data.

    Only generates the angular portion of the sphere near the actual
    nuclei, so a thin SPIM cap doesn't produce a giant full sphere.

    Parameters
    ----------
    center : (3,) sphere center (x, y, z)
    radius : sphere radius
    data_points : (N, 3) positions of nuclei (x, y, z)
    padding_deg : extra angular margin around the data
    """
    d = data_points - center
    r_safe = np.maximum(np.sqrt(np.sum(d**2, axis=1)), 1e-10)

    # Spherical: theta from -Y axis (AP at top in napari), phi in XZ
    theta = np.arccos(np.clip(-d[:, 1] / r_safe, -1, 1))
    phi = np.arctan2(d[:, 2], d[:, 0])

    pad = np.radians(padding_deg)
    theta_min = max(0, np.percentile(theta, 2) - pad)
    theta_max = min(np.pi, np.percentile(theta, 98) + pad)
    phi_min = np.percentile(phi, 2) - pad
    phi_max = np.percentile(phi, 98) + pad

    tg, pg = np.meshgrid(
        np.linspace(theta_min, theta_max, n_lat),
        np.linspace(phi_min, phi_max, n_lon),
        indexing="ij",
    )

    x = center[0] + radius * np.sin(tg) * np.cos(pg)
    y = center[1] - radius * np.cos(tg)  # AP at -Y (top in napari)
    z = center[2] + radius * np.sin(tg) * np.sin(pg)

    verts = np.column_stack([x.ravel(), y.ravel(), z.ravel()])

    faces = []
    for i in range(n_lat - 1):
        for j in range(n_lon - 1):
            v0 = i * n_lon + j
            v1 = v0 + 1
            v2 = v0 + n_lon
            v3 = v2 + 1
            faces.append([v0, v1, v3])
            faces.append([v0, v3, v2])

    return verts, np.array(faces)


# ############################################################################
#                     INGRESSION ANALYSIS
# ############################################################################


def compute_radial_velocity(
    df: pd.DataFrame, center: np.ndarray, smooth_window: int = 5
) -> pd.DataFrame:
    """Compute per-spot radial velocity with smoothing and enhanced track stats.

    Negative = moving inward (ingressing).
    Positive = moving outward.

    Parameters
    ----------
    df : spots DataFrame with POSITION_X/Y/Z, TRACK_ID, FRAME
    center : (3,) sphere center
    smooth_window : rolling-window size for smoothed V_r (odd, ≥3)

    Returns
    -------
    DataFrame with new columns:
      RADIAL_DIST_TO_CENTER   – distance from sphere centre
      RADIAL_VELOCITY         – instantaneous V_r (µm/frame)
      RADIAL_VELOCITY_SMOOTH  – smoothed V_r (rolling mean per track)
      TRACK_MEDIAN_RADIAL_VEL – per-track median of smoothed V_r (robust)
      TRACK_MEAN_RADIAL_VEL   – per-track mean of smoothed V_r
      TRACK_INWARD_FRACTION   – fraction of track points with V_r < 0
      TRACK_MAX_INWARD_VEL    – most negative V_r per track (peak inward)
    """
    x = df["POSITION_X"].values - center[0]
    y = df["POSITION_Y"].values - center[1]
    z = df["POSITION_Z"].values - center[2]
    rdist = np.sqrt(x**2 + y**2 + z**2)
    df["RADIAL_DIST_TO_CENTER"] = rdist

    nan_cols = {
        "RADIAL_VELOCITY": np.nan,
        "RADIAL_VELOCITY_SMOOTH": np.nan,
        "TRACK_MEDIAN_RADIAL_VEL": np.nan,
        "TRACK_MEAN_RADIAL_VEL": np.nan,
        "TRACK_INWARD_FRACTION": np.nan,
        "TRACK_MAX_INWARD_VEL": np.nan,
    }
    if "TRACK_ID" not in df.columns:
        for col in nan_cols:
            df[col] = np.full(len(df), np.nan)
        return df

    # --- Instantaneous radial velocity ---
    sorted_idx = df.sort_values(["TRACK_ID", "FRAME"]).index
    tid = df.loc[sorted_idx, "TRACK_ID"].values
    r = rdist[sorted_idx]
    f = df.loc[sorted_idx, "FRAME"].values.astype(float)

    dr = np.diff(r, prepend=np.nan)
    dt = np.diff(f, prepend=-9999)
    same_track = np.diff(tid.astype(float), prepend=-9999) == 0
    vel = np.where(same_track & (dt > 0), dr / dt, np.nan)

    rv = np.full(len(df), np.nan)
    rv[sorted_idx] = vel
    df["RADIAL_VELOCITY"] = rv

    # --- Smoothed V_r: rolling mean within each track ---
    sw = max(3, smooth_window | 1)  # ensure odd ≥3
    tmp = pd.DataFrame(
        {
            "TRACK_ID": df["TRACK_ID"],
            "FRAME": df["FRAME"],
            "RV": rv,
            "_idx": df.index,
        }
    ).sort_values(["TRACK_ID", "FRAME"])

    tmp["RV_SMOOTH"] = tmp.groupby("TRACK_ID")["RV"].transform(
        lambda s: s.rolling(sw, center=True, min_periods=1).mean()
    )
    rv_smooth = np.full(len(df), np.nan)
    rv_smooth[tmp["_idx"].values] = tmp["RV_SMOOTH"].values
    df["RADIAL_VELOCITY_SMOOTH"] = rv_smooth

    # --- Per-track summary statistics (based on smoothed V_r) ---
    valid_mask = ~np.isnan(rv_smooth)
    if valid_mask.any():
        sub = df.loc[valid_mask, ["TRACK_ID", "RADIAL_VELOCITY_SMOOTH"]]
        grp = sub.groupby("TRACK_ID")["RADIAL_VELOCITY_SMOOTH"]
        med_map = grp.median()
        mean_map = grp.mean()
        inward_frac_map = grp.agg(lambda s: (s < 0).mean())
        max_inward_map = grp.min()  # most negative = strongest inward

        df["TRACK_MEDIAN_RADIAL_VEL"] = df["TRACK_ID"].map(med_map)
        df["TRACK_MEAN_RADIAL_VEL"] = df["TRACK_ID"].map(mean_map)
        df["TRACK_INWARD_FRACTION"] = df["TRACK_ID"].map(inward_frac_map)
        df["TRACK_MAX_INWARD_VEL"] = df["TRACK_ID"].map(max_inward_map)
    else:
        for col in [
            "TRACK_MEDIAN_RADIAL_VEL",
            "TRACK_MEAN_RADIAL_VEL",
            "TRACK_INWARD_FRACTION",
            "TRACK_MAX_INWARD_VEL",
        ]:
            df[col] = np.nan

    return df


def flag_ingressing(
    df: pd.DataFrame,
    threshold: float = -0.5,
    min_track_frames: int = 5,
    min_inward_fraction: float = 0.5,
) -> pd.DataFrame:
    """Flag tracks as ingressing using robust multi-criteria detection.

    A track is flagged as INGRESSING if ALL of:
      1. Track has ≥ *min_track_frames* points
      2. Track **median** smoothed V_r ≤ *threshold*
      3. ≥ *min_inward_fraction* of track points have V_r < 0

    Parameters
    ----------
    df : spots DataFrame with TRACK_MEDIAN_RADIAL_VEL, TRACK_INWARD_FRACTION
    threshold : µm/frame – tracks with median V_r ≤ threshold → ingressing
    min_track_frames : minimum track length
    min_inward_fraction : minimum fraction of inward-moving points

    Returns
    -------
    DataFrame with INGRESSING (bool) column
    """
    has_med = "TRACK_MEDIAN_RADIAL_VEL" in df.columns
    has_mean = "TRACK_MEAN_RADIAL_VEL" in df.columns
    if "TRACK_ID" not in df.columns or (not has_med and not has_mean):
        df["INGRESSING"] = False
        return df

    # Use median if available (more robust), else fall back to mean
    vel_col = "TRACK_MEDIAN_RADIAL_VEL" if has_med else "TRACK_MEAN_RADIAL_VEL"
    vel = df[vel_col].values

    track_len = df.groupby("TRACK_ID")["FRAME"].transform("count")
    long_enough = track_len >= min_track_frames

    # Inward fraction criterion
    if "TRACK_INWARD_FRACTION" in df.columns:
        inward_ok = df["TRACK_INWARD_FRACTION"].values >= min_inward_fraction
    else:
        inward_ok = np.ones(len(df), dtype=bool)

    df["INGRESSING"] = long_enough & (~np.isnan(vel)) & (vel <= threshold) & inward_ok
    return df


def auto_ingression_threshold(df: pd.DataFrame) -> tuple[float, str]:
    """Auto-detect ingression V_r threshold using Otsu's method.

    Finds the threshold on the distribution of per-track median V_r
    that maximises between-class variance (ingressing vs non-ingressing).

    Returns
    -------
    (threshold, description_string)
    """
    col = "TRACK_MEDIAN_RADIAL_VEL"
    if col not in df.columns or "TRACK_ID" not in df.columns:
        return -0.5, "default (no data)"

    # Unique per-track values
    track_stats = df.dropna(subset=[col]).groupby("TRACK_ID")[col].first()
    vals = track_stats.values
    vals = vals[~np.isnan(vals)]

    if len(vals) < 20:
        return float(np.median(vals)) if len(vals) > 0 else -0.5, "median (few tracks)"

    # Otsu 1D: maximise between-class variance
    n_bins = 200
    hist, edges = np.histogram(vals, bins=n_bins)
    centres = (edges[:-1] + edges[1:]) / 2
    total = hist.sum()
    if total == 0:
        return -0.5, "default (empty)"

    best_thresh, best_var = centres[0], 0.0
    w0, sum0 = 0, 0.0
    sum_total = float((hist * centres).sum())

    for i in range(n_bins):
        w0 += hist[i]
        if w0 == 0:
            continue
        w1 = total - w0
        if w1 == 0:
            break
        sum0 += hist[i] * centres[i]
        m0 = sum0 / w0
        m1 = (sum_total - sum0) / w1
        var = w0 * w1 * (m0 - m1) ** 2
        if var > best_var:
            best_var = var
            best_thresh = float(centres[i])

    # Sanity: for ingression we expect negative thresholds
    best_thresh = min(best_thresh, -0.005)
    p_below = float((vals <= best_thresh).mean() * 100)

    desc = f"Otsu={best_thresh:.3f} µm/f ({p_below:.0f}% of tracks below)"
    return best_thresh, desc


# ############################################################################
#                        TRACKS → napari format
# ############################################################################


def spots_to_napari_tracks(
    df: pd.DataFrame,
    max_tracks: int = 5000,
    min_length: int = 3,
    frame_min: int | None = None,
    frame_max: int | None = None,
) -> np.ndarray | None:
    """Convert spots DataFrame to napari Tracks format.

    napari Tracks format: [track_id, t, z, y, x] — 5 columns for 3D+T.

    Parameters
    ----------
    frame_min, frame_max : optional frame range filter for tracks.
        Only track points within [frame_min, frame_max] are included.
        Tracks that have no points in range are dropped entirely.
    """
    if "TRACK_ID" not in df.columns:
        return None

    df_t = df.dropna(subset=["TRACK_ID"])
    if df_t.empty:
        return None

    # Frame range filter
    if frame_min is not None:
        df_t = df_t[df_t["FRAME"] >= frame_min]
    if frame_max is not None:
        df_t = df_t[df_t["FRAME"] <= frame_max]
    if df_t.empty:
        return None

    # Filter by min length (within the frame window)
    track_len = df_t.groupby("TRACK_ID")["FRAME"].transform("count")
    df_t = df_t[track_len >= min_length]

    # Subsample tracks if needed
    tids = df_t["TRACK_ID"].unique()
    if len(tids) > max_tracks:
        rng = np.random.default_rng(42)
        tids = rng.choice(tids, max_tracks, replace=False)
        df_t = df_t[df_t["TRACK_ID"].isin(tids)]

    df_t = df_t.sort_values(["TRACK_ID", "FRAME"])

    # napari tracks: [ID, t, z, y, x]
    tracks = np.column_stack(
        [
            df_t["TRACK_ID"].values,
            df_t["FRAME"].values,
            df_t["POSITION_Z"].values,
            df_t["POSITION_Y"].values,
            df_t["POSITION_X"].values,
        ]
    )
    return tracks


# ############################################################################
#                            NAPARI VIEWER
# ############################################################################


class EmbryoViewer:
    """napari-based 4D embryo viewer with orientation & sphere fitting."""

    def __init__(
        self,
        spots_df: pd.DataFrame | None = None,
        tracks_df: pd.DataFrame | None = None,
    ):
        napari, magicgui = _import_napari()

        self.spots = spots_df
        self.tracks_raw = tracks_df

        # State
        self.animal_pole: np.ndarray | None = None
        self.dorsal_mark: np.ndarray | None = None
        self.oriented = False
        self.sphere_fit_result: dict | None = None
        self.transform_params: dict | None = None
        self.spots_layer = None
        self.tracks_layer = None
        self._tracks_layer_secondary = None  # for dual-layer ingression view
        self._display_pct: float = 100.0
        self._max_tracks: int = 5000
        self._track_frame_min: int | None = None
        self._track_frame_max: int | None = None
        self._roi_bounds: dict | None = None  # {x_min, x_max, y_min, ...}
        self._roi_only: bool = False  # show only ROI spots
        self._roi_shapes_layer = None
        self._roi_updating: bool = False  # guard for event recursion
        self._drawing_roi: bool = False  # guard: suppress rebuilds while drawing
        self._current_color_mode: str = "track_id"  # tracks colour mode
        self._tracks_cache_key: tuple | None = None  # for caching track builds

        # Create viewer
        self.viewer = napari.Viewer(title="Embryo 4D Viewer", ndisplay=3)

        # Load data layers if provided
        if self.spots is not None:
            self._add_spots_layer()
            self._add_tracks_layer()

        # ── Landmarks layer (3D so visible at every frame) ──
        self.landmarks_layer = self.viewer.add_points(
            np.empty((0, 3)),
            name="Landmarks",
            size=20,
            face_color="red",
            symbol="diamond",
            ndim=3,
            border_width=0,
        )

        # ── Build control widgets ──
        self._build_widgets()

        # Connect click handler for landmark selection
        self.viewer.mouse_drag_callbacks.append(self._on_click)

    def _add_spots_layer(self):
        """Add the main nuclei points layer."""
        data, colors = self._build_display_points()
        self.spots_layer = self.viewer.add_points(
            data,
            name="Nuclei",
            size=2,
            face_color=colors,
            border_width=0,
            ndim=4,
        )

    def _build_display_points(self):
        """Build display arrays, respecting display %, ROI and ingression filters."""
        df = self.spots
        n = len(df)
        mask = np.ones(n, dtype=bool)

        # ROI filter
        if self._roi_only and self._roi_bounds is not None:
            b = self._roi_bounds
            mask &= (
                (df["POSITION_X"].values >= b["x_min"])
                & (df["POSITION_X"].values <= b["x_max"])
                & (df["POSITION_Y"].values >= b["y_min"])
                & (df["POSITION_Y"].values <= b["y_max"])
            )

        # Ingression filter (when "show ingressing tracks only" is checked)
        _ing_w = getattr(self, "ingression_tracks_only", None)
        if _ing_w is not None and _ing_w.value and "INGRESSING" in df.columns:
            ing_tids = df.loc[df["INGRESSING"].astype(bool), "TRACK_ID"].unique()
            mask &= df["TRACK_ID"].isin(ing_tids).values

        pool = np.where(mask)[0]

        pct = self._display_pct
        if pct < 100:
            n_show = max(1000, int(len(pool) * pct / 100))
            if n_show < len(pool):
                rng = np.random.default_rng(0)
                idx = rng.choice(pool, n_show, replace=False)
                idx.sort()
            else:
                idx = pool
        else:
            idx = pool

        data = np.column_stack(
            [
                df["FRAME"].values[idx],
                df["POSITION_Z"].values[idx],
                df["POSITION_Y"].values[idx],
                df["POSITION_X"].values[idx],
            ]
        )
        frame = df["FRAME"].values[idx].astype(float)
        frame_norm = (frame - frame.min()) / max(np.ptp(frame), 1)
        return data, _viridis_colors(frame_norm)

    def _add_tracks_layer(self):
        """Add tracks as a napari Tracks layer."""
        max_t = int(self._max_tracks)
        tracks = spots_to_napari_tracks(
            self.spots,
            max_tracks=max_t,
            frame_min=self._track_frame_min,
            frame_max=self._track_frame_max,
        )
        if tracks is not None and len(tracks) > 0:
            self.tracks_layer = self.viewer.add_tracks(
                tracks,
                name="Tracks",
                tail_width=1,
                tail_length=40,
                color_by="track_id",
            )
        else:
            self.tracks_layer = None

    def _build_tracks_data(self):
        """Build tracks array + filtered DataFrame applying all active filters.

        Respects: ROI-only, ingression-only, time window, max tracks.
        Returns (tracks_array, df_filtered) or (None, None).
        """
        if self.spots is None or "TRACK_ID" not in self.spots.columns:
            return None, None

        df = self.spots.dropna(subset=["TRACK_ID"]).copy()
        if df.empty:
            return None, None

        # ROI filter: keep tracks that have at least one spot in ROI
        if self._roi_only and "IN_ROI" in df.columns:
            roi_tids = df.loc[df["IN_ROI"].astype(bool), "TRACK_ID"].unique()
            df = df[df["TRACK_ID"].isin(roi_tids)]

        # Ingression filter
        _ing_w = getattr(self, "ingression_tracks_only", None)
        if _ing_w is not None and _ing_w.value and "INGRESSING" in df.columns:
            ing_tids = df.loc[df["INGRESSING"].astype(bool), "TRACK_ID"].unique()
            df = df[df["TRACK_ID"].isin(ing_tids)]

        # Time window
        if self._track_frame_min is not None:
            df = df[df["FRAME"] >= self._track_frame_min]
        if self._track_frame_max is not None:
            df = df[df["FRAME"] <= self._track_frame_max]
        if df.empty:
            return None, None

        # Min track length
        track_len = df.groupby("TRACK_ID")["FRAME"].transform("count")
        df = df[track_len >= 3]
        if df.empty:
            return None, None

        # Subsample tracks
        tids = df["TRACK_ID"].unique()
        max_t = int(self._max_tracks)
        if len(tids) > max_t:
            rng = np.random.default_rng(42)
            tids = rng.choice(tids, max_t, replace=False)
            df = df[df["TRACK_ID"].isin(tids)]

        df = df.sort_values(["TRACK_ID", "FRAME"])
        tracks = np.column_stack(
            [
                df["TRACK_ID"].values,
                df["FRAME"].values,
                df["POSITION_Z"].values,
                df["POSITION_Y"].values,
                df["POSITION_X"].values,
            ]
        )
        return tracks, df

    def _rebuild_tracks(self):
        """Rebuild the tracks layer(s) with all active filters + current colour.

        In 'ingressing' colour mode, creates TWO layers:
          - "Other Tracks" (grey, thin, translucent)
          - "Ingressing Tracks" (magenta, thick, opaque)
        Otherwise creates a single "Tracks" layer.
        """
        if self._drawing_roi:
            return  # suppress rebuilds during ROI drawing
        self._remove_layer("tracks_layer")
        self._remove_layer("_tracks_layer_secondary")

        result = self._build_tracks_data()
        if result is None or result[0] is None:
            return
        tracks, df_t = result
        if len(tracks) == 0:
            return

        mode = self._current_color_mode

        # ── Dual-layer mode for ingression classification view ──
        if mode == "ingressing" and "INGRESSING" in df_t.columns:
            ing_mask = df_t["INGRESSING"].values.astype(bool)
            # Non-ingressing tracks
            if (~ing_mask).any():
                df_other = df_t[~ing_mask]
                if len(df_other) >= 3:
                    # filter to tracks with ≥2 points
                    tc = df_other.groupby("TRACK_ID")["FRAME"].transform("count")
                    df_other = df_other[tc >= 2].sort_values(["TRACK_ID", "FRAME"])
                    if len(df_other) > 0:
                        t_other = np.column_stack(
                            [
                                df_other["TRACK_ID"].values,
                                df_other["FRAME"].values,
                                df_other["POSITION_Z"].values,
                                df_other["POSITION_Y"].values,
                                df_other["POSITION_X"].values,
                            ]
                        )
                        self._tracks_layer_secondary = self.viewer.add_tracks(
                            t_other,
                            name="Other Tracks",
                            tail_width=0.5,
                            tail_length=30,
                            color_by="track_id",
                            colormap="gray",
                            opacity=0.25,
                        )
            # Ingressing tracks
            if ing_mask.any():
                df_ing = df_t[ing_mask]
                tc = df_ing.groupby("TRACK_ID")["FRAME"].transform("count")
                df_ing = df_ing[tc >= 2].sort_values(["TRACK_ID", "FRAME"])
                if len(df_ing) > 0:
                    t_ing = np.column_stack(
                        [
                            df_ing["TRACK_ID"].values,
                            df_ing["FRAME"].values,
                            df_ing["POSITION_Z"].values,
                            df_ing["POSITION_Y"].values,
                            df_ing["POSITION_X"].values,
                        ]
                    )
                    self.tracks_layer = self.viewer.add_tracks(
                        t_ing,
                        name="Ingressing Tracks",
                        tail_width=2,
                        tail_length=50,
                        color_by="track_id",
                        colormap="magenta",
                    )
            return

        # ── Single-layer mode (all other colour modes) ──
        col_map = {
            "depth": "SPHERICAL_DEPTH",
            "theta": "THETA_DEG",
            "phi": "PHI_DEG",
            "radial_vel": "RADIAL_VELOCITY_SMOOTH",
            "frame": "FRAME",
        }
        props = {}
        color_by = "track_id"
        cmap = "hsv"

        if mode == "speed":
            sp = self._compute_track_speed(df_t)
            nan_m = np.isnan(sp)
            if not nan_m.all():
                sp[nan_m] = 0.0
                props["color_value"] = sp
                color_by = "color_value"
                cmap = "viridis"
        elif mode in col_map and col_map[mode] in df_t.columns:
            vals = df_t[col_map[mode]].values.astype(float)
            nan_m = np.isnan(vals)
            if not nan_m.all():
                vals[nan_m] = np.nanmedian(vals)
                props["color_value"] = vals
                color_by = "color_value"
                cmap = "turbo" if mode == "radial_vel" else "viridis"

        kw = dict(
            name="Tracks",
            tail_width=1,
            tail_length=40,
            color_by=color_by,
            colormap=cmap,
        )
        if props:
            kw["properties"] = props
        self.tracks_layer = self.viewer.add_tracks(tracks, **kw)

    @staticmethod
    def _compute_track_speed(df_t: pd.DataFrame) -> np.ndarray:
        """Per-point speed (\u00b5m/frame) for a sorted tracks DataFrame."""
        tid = df_t["TRACK_ID"].values
        x, y, z = (
            df_t["POSITION_X"].values,
            df_t["POSITION_Y"].values,
            df_t["POSITION_Z"].values,
        )
        f = df_t["FRAME"].values.astype(float)
        dx = np.diff(x, prepend=np.nan)
        dy = np.diff(y, prepend=np.nan)
        dz = np.diff(z, prepend=np.nan)
        df_v = np.diff(f, prepend=-9999)
        same = np.diff(tid.astype(float), prepend=-9999) == 0
        dist = np.sqrt(dx**2 + dy**2 + dz**2)
        dt = np.maximum(df_v, 1.0)
        return np.where(same & (df_v > 0), dist / dt, np.nan)

    def _remove_layer(self, attr_name: str):
        """Safely remove a napari layer stored as self.<attr_name>."""
        layer = getattr(self, attr_name, None)
        if layer is not None:
            try:
                self.viewer.layers.remove(layer)
            except (ValueError, KeyError):
                pass
            setattr(self, attr_name, None)

    def _reset_downstream(self):
        """Clear orientation / sphere / ROI state after new data load."""
        self.animal_pole = None
        self.dorsal_mark = None
        self.oriented = False
        self.sphere_fit_result = None
        self.transform_params = None
        self._roi_bounds = None
        self._roi_only = False
        self.roi_show_only.value = False
        if self._roi_shapes_layer is not None:
            try:
                self.viewer.layers.remove(self._roi_shapes_layer)
            except (ValueError, KeyError):
                pass
            self._roi_shapes_layer = None
        self.lbl_ap.value = "AP: not set"
        self.lbl_dorsal.value = "Dorsal: not set"
        self.lbl_sphere.value = "Sphere: not fitted"
        self.lbl_roi.value = "ROI: not set"
        self.lbl_ingression.value = "Ingression: not computed"
        self.ingression_tracks_only.value = False

        # Update frame range sliders to new data
        if self.spots is not None:
            fmin = int(self.spots["FRAME"].min())
            fmax = int(self.spots["FRAME"].max())
            self.track_frame_start.min = fmin
            self.track_frame_start.max = fmax
            self.track_frame_start.value = fmin
            self.track_frame_end.min = fmin
            self.track_frame_end.max = fmax
            self.track_frame_end.value = fmax
            self._track_frame_min = None
            self._track_frame_max = None

    def _load_spots_from_ui(self):
        """Open a file dialog to load a spots CSV (background thread)."""
        from qtpy.QtWidgets import QFileDialog
        from qtpy.QtCore import QTimer

        spots_path, _ = QFileDialog.getOpenFileName(
            self.viewer.window._qt_window,
            "Select Spots CSV",
            "",
            "CSV files (*.csv);;All files (*)",
        )
        if not spots_path:
            return

        self.lbl_spots_status.value = "\u23f3 Loading (please wait)..."
        self.viewer.status = "Loading spots in background thread..."

        self._pending_spots_df = None
        self._pending_spots_err = None

        def _bg():
            try:
                self._pending_spots_df = load_spots(spots_path)
            except Exception as e:
                self._pending_spots_err = str(e)

        t = threading.Thread(target=_bg, daemon=True)
        t.start()

        def _poll():
            if t.is_alive():
                QTimer.singleShot(250, _poll)
                return
            if self._pending_spots_err:
                self.lbl_spots_status.value = f"Error: {self._pending_spots_err}"
                self.viewer.status = f"Load failed: {self._pending_spots_err}"
                return
            self._finish_spots_load(self._pending_spots_df)

        QTimer.singleShot(250, _poll)

    def _finish_spots_load(self, df):
        """Finalise spots load on the main thread."""
        self.spots = df
        self._remove_layer("tracks_layer")
        self._remove_layer("_tracks_layer_secondary")
        self._remove_layer("spots_layer")

        self.viewer.status = "Building point cloud..."
        self._add_spots_layer()
        self._add_tracks_layer()

        n_s = f"{len(self.spots):,}"
        n_f = self.spots["FRAME"].nunique()
        has_tid = (
            "TRACK_ID" in self.spots.columns and self.spots["TRACK_ID"].notna().any()
        )
        msg = f"{n_s} spots, {n_f} frames"
        if self.tracks_layer is not None:
            msg += " + trajectories"
        elif not has_tid:
            msg += " (no TRACK_ID)"
        self.lbl_spots_status.value = msg
        self.viewer.status = msg + " | Use time slider at bottom to scrub frames"
        self._reset_downstream()
        self.viewer.reset_view()

    def _load_tracks_from_ui(self):
        """Open a file dialog to load a tracks summary CSV."""
        from qtpy.QtWidgets import QFileDialog, QApplication

        tracks_path, _ = QFileDialog.getOpenFileName(
            self.viewer.window._qt_window,
            "Select Tracks CSV",
            "",
            "CSV files (*.csv);;All files (*)",
        )
        if not tracks_path:
            return

        self.viewer.status = "Loading tracks summary..."
        QApplication.processEvents()
        try:
            self.tracks_raw = load_trackmate_csv(tracks_path)
        except Exception as e:
            self.viewer.status = f"Error loading tracks: {e}"
            self.lbl_tracks_status.value = f"Error: {e}"
            return

        n_t = f"{len(self.tracks_raw):,}"
        cols = [c for c in self.tracks_raw.columns if "TRACK" in c.upper()]
        msg = f"{n_t} track records"
        self.lbl_tracks_status.value = msg
        self.viewer.status = msg

    def _build_widgets(self):
        """Create magicgui dock widgets for controls."""
        napari, magicgui = _import_napari()
        from magicgui import magicgui as mg_deco
        from magicgui.widgets import (
            PushButton,
            Label,
            FloatSlider,
            CheckBox,
            Container,
        )

        # Determine frame range for sliders
        if self.spots is not None:
            fmin_data = int(self.spots["FRAME"].min())
            fmax_data = int(self.spots["FRAME"].max())
        else:
            fmin_data, fmax_data = 0, 1000

        # ── File loading (separate buttons) ──
        btn_load_spots = PushButton(text="Load Spots CSV")
        btn_load_spots.changed.connect(self._load_spots_from_ui)
        self.lbl_spots_status = Label(
            value=f"{len(self.spots):,} spots loaded"
            if self.spots is not None
            else "No spots loaded"
        )

        btn_load_tracks = PushButton(text="Load Tracks CSV")
        btn_load_tracks.changed.connect(self._load_tracks_from_ui)
        self.lbl_tracks_status = Label(
            value=f"{len(self.tracks_raw):,} track records"
            if self.tracks_raw is not None
            else "No tracks loaded"
        )

        # ── Display controls ──
        self.display_pct_slider = FloatSlider(
            value=100, min=1, max=100, step=1, label="Display %"
        )
        self.display_pct_slider.changed.connect(self._on_display_pct_changed)

        self.max_tracks_slider = FloatSlider(
            value=5000, min=100, max=20000, step=500, label="Max tracks"
        )
        self.max_tracks_slider.changed.connect(self._on_max_tracks_changed)

        # ── Track time range ──
        self.track_frame_start = FloatSlider(
            value=fmin_data, min=fmin_data, max=fmax_data, step=1, label="▶ Start frame"
        )
        self.track_frame_end = FloatSlider(
            value=fmax_data, min=fmin_data, max=fmax_data, step=1, label="■ End frame"
        )
        self.track_frame_start.changed.connect(self._on_track_frame_changed)
        self.track_frame_end.changed.connect(self._on_track_frame_changed)

        # ── Landmark picking ──
        self._pick_mode = "none"

        btn_ap = PushButton(text="Pick Animal Pole")
        btn_ap.changed.connect(lambda: self._set_pick_mode("ap"))

        btn_dorsal = PushButton(text="Pick Dorsal")
        btn_dorsal.changed.connect(lambda: self._set_pick_mode("dorsal"))

        self.lbl_ap = Label(value="AP: not set")
        self.lbl_dorsal = Label(value="Dorsal: not set")

        # ── Orient ──
        btn_orient = PushButton(text="Orient Embryo")
        btn_orient.changed.connect(self._orient)

        # ── Sphere fit ──
        btn_sphere = PushButton(text="Fit Sphere")
        btn_sphere.changed.connect(self._fit_sphere)
        self.lbl_sphere = Label(value="Sphere: not fitted")

        # ── ROI rectangle selection ──
        btn_draw_roi = PushButton(text="▣ Draw ROI Rectangle")
        btn_draw_roi.changed.connect(self._activate_roi_drawing)
        self.lbl_roi = Label(value="ROI: not set")

        self.roi_show_only = CheckBox(value=False, text="Show ROI spots only")
        self.roi_show_only.changed.connect(self._on_roi_toggle)

        btn_clear_roi = PushButton(text="Clear ROI")
        btn_clear_roi.changed.connect(self._clear_roi)

        # ── Ingression analysis ──
        btn_compute_ingression = PushButton(text="\u2b07 Compute Radial Velocity")
        btn_compute_ingression.changed.connect(self._compute_radial_velocity)
        self.ingression_threshold = FloatSlider(
            value=-0.1,
            min=-2.0,
            max=0.5,
            step=0.01,
            label="V_r threshold",
        )
        self.ingression_inward_frac = FloatSlider(
            value=0.5,
            min=0.0,
            max=1.0,
            step=0.05,
            label="Min inward frac",
        )
        self.ingression_min_frames = FloatSlider(
            value=5,
            min=2,
            max=30,
            step=1,
            label="Min track frames",
        )
        btn_auto_threshold = PushButton(text="\U0001f3af Auto-Detect Threshold")
        btn_auto_threshold.changed.connect(self._auto_detect_threshold)
        btn_flag_ingression = PushButton(text="Flag Ingressing Cells")
        btn_flag_ingression.changed.connect(self._flag_ingression)
        self.lbl_ingression = Label(value="Ingression: not computed")

        self.ingression_tracks_only = CheckBox(
            value=False, text="Show ingressing tracks only"
        )
        self.ingression_tracks_only.changed.connect(self._on_ingression_tracks_toggle)

        # ── Colour mode ──
        btn_color_frame = PushButton(text="Colour by Frame")
        btn_color_frame.changed.connect(lambda: self._recolor("frame"))
        btn_color_depth = PushButton(text="Colour by Depth")
        btn_color_depth.changed.connect(lambda: self._recolor("depth"))
        btn_color_theta = PushButton(text="Colour by Latitude")
        btn_color_theta.changed.connect(lambda: self._recolor("theta"))
        btn_color_phi = PushButton(text="Colour by Longitude")
        btn_color_phi.changed.connect(lambda: self._recolor("phi"))
        btn_color_speed = PushButton(text="Colour by Speed")
        btn_color_speed.changed.connect(lambda: self._recolor("speed"))
        btn_color_radvel = PushButton(text="Colour by Radial Velocity")
        btn_color_radvel.changed.connect(lambda: self._recolor("radial_vel"))
        btn_color_ingress = PushButton(text="Colour Ingressing")
        btn_color_ingress.changed.connect(lambda: self._recolor("ingressing"))
        btn_color_roi = PushButton(text="Colour ROI")
        btn_color_roi.changed.connect(lambda: self._recolor("roi"))

        # ── Export & Video ──
        btn_export = PushButton(text="Export Enriched CSV")
        btn_export.changed.connect(self._export)
        self.lbl_export = Label(value="")
        btn_record_video = PushButton(text="🎬 Record Time-lapse Video")
        btn_record_video.changed.connect(self._record_video)

        container = Container(
            widgets=[
                Label(value="── Data ──"),
                btn_load_spots,
                self.lbl_spots_status,
                btn_load_tracks,
                self.lbl_tracks_status,
                Label(value="━━ ⏱ TIME WINDOW ━━"),
                self.track_frame_start,
                self.track_frame_end,
                Label(value="── Display ──"),
                self.display_pct_slider,
                self.max_tracks_slider,
                Label(value="── Landmarks ──"),
                btn_ap,
                self.lbl_ap,
                btn_dorsal,
                self.lbl_dorsal,
                Label(value="── Orientation ──"),
                btn_orient,
                Label(value="── Sphere Fit ──"),
                btn_sphere,
                self.lbl_sphere,
                Label(value="── ROI Selection ──"),
                btn_draw_roi,
                self.lbl_roi,
                self.roi_show_only,
                btn_clear_roi,
                Label(value="━━ ⬇ INGRESSION ━━"),
                btn_compute_ingression,
                self.ingression_threshold,
                self.ingression_inward_frac,
                self.ingression_min_frames,
                btn_auto_threshold,
                btn_flag_ingression,
                self.lbl_ingression,
                self.ingression_tracks_only,
                Label(value="── Colouring ──"),
                btn_color_frame,
                btn_color_depth,
                btn_color_theta,
                btn_color_phi,
                btn_color_speed,
                btn_color_radvel,
                btn_color_ingress,
                btn_color_roi,
                Label(value="── Export & Video ──"),
                btn_export,
                btn_record_video,
                self.lbl_export,
            ]
        )

        # Wrap in a QScrollArea so the panel is scrollable on small screens
        from qtpy.QtWidgets import QScrollArea, QWidget, QVBoxLayout
        from qtpy.QtCore import Qt

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setMinimumWidth(300)
        inner = container.native  # the underlying QWidget from magicgui
        scroll.setWidget(inner)

        self.viewer.window.add_dock_widget(scroll, name="Controls", area="right")

    # ── Landmark picking ──────────────────────────────────────────────

    def _set_pick_mode(self, mode: str):
        if self._pick_mode == mode:
            self._pick_mode = "none"
            self.viewer.status = "Pick mode OFF"
        else:
            self._pick_mode = mode
            labels = {
                "ap": "ANIMAL POLE",
                "dorsal": "DORSAL",
            }
            label = labels.get(mode, mode.upper())
            self.viewer.status = (
                f"\U0001f3af CLICK near the {label} \u2014 "
                f"snaps to nearest spot in current frame. "
                f"Click button again to cancel."
            )

    def _on_click(self, viewer, event):
        """Handle mouse click for landmark picking (snaps to nearest spot)."""
        if self._pick_mode == "none":
            return
        if self.spots_layer is None or self.spots is None:
            self.viewer.status = "Load spots first!"
            return

        # Get click world coords \u2014 napari gives (t, z, y, x) for 4D
        world = np.array(event.position)
        if len(world) < 4:
            return
        click_x, click_y, click_z = world[3], world[2], world[1]

        # Find nearest spot in the current frame
        current_frame = int(round(self.viewer.dims.current_step[0]))
        frame_mask = self.spots["FRAME"].values == current_frame
        if not frame_mask.any():
            frames = self.spots["FRAME"].unique()
            current_frame = int(frames[np.argmin(np.abs(frames - current_frame))])
            frame_mask = self.spots["FRAME"].values == current_frame

        sx = self.spots.loc[frame_mask, "POSITION_X"].values
        sy = self.spots.loc[frame_mask, "POSITION_Y"].values
        sz = self.spots.loc[frame_mask, "POSITION_Z"].values
        dists = np.sqrt((sx - click_x) ** 2 + (sy - click_y) ** 2 + (sz - click_z) ** 2)
        nearest = int(dists.argmin())
        pos_xyz = np.array([sx[nearest], sy[nearest], sz[nearest]])
        snap_d = dists[nearest]

        if self._pick_mode == "ap":
            self.animal_pole = pos_xyz
            self.lbl_ap.value = (
                f"AP: ({pos_xyz[0]:.1f}, {pos_xyz[1]:.1f}, {pos_xyz[2]:.1f})"
            )
            self._pick_mode = "none"
            self.viewer.status = (
                f"Animal Pole set (snapped {snap_d:.1f}\u00b5m to nearest spot)"
            )
        elif self._pick_mode == "dorsal":
            self.dorsal_mark = pos_xyz
            self.lbl_dorsal.value = (
                f"Dorsal: ({pos_xyz[0]:.1f}, {pos_xyz[1]:.1f}, {pos_xyz[2]:.1f})"
            )
            self._pick_mode = "none"
            self.viewer.status = (
                f"Dorsal set (snapped {snap_d:.1f}\u00b5m to nearest spot)"
            )
        self._update_landmarks()

    def _update_landmarks(self):
        """Update landmarks (3D layer \u2014 visible at every frame)."""
        pts = []
        colors = []
        if self.animal_pole is not None:
            # 3D: (z, y, x) \u2014 visible at all time steps
            pts.append([self.animal_pole[2], self.animal_pole[1], self.animal_pole[0]])
            colors.append("red")
        if self.dorsal_mark is not None:
            pts.append([self.dorsal_mark[2], self.dorsal_mark[1], self.dorsal_mark[0]])
            colors.append("cyan")

        if pts:
            self.landmarks_layer.data = np.array(pts)
            self.landmarks_layer.face_color = colors
            self.landmarks_layer.size = 20
        else:
            self.landmarks_layer.data = np.empty((0, 3))

    # ── Orientation ───────────────────────────────────────────────────

    def _orient(self):
        """Apply AP → +Y, Dorsal → +X rotation."""
        if self.animal_pole is None or self.dorsal_mark is None:
            self.viewer.status = "Set both Animal Pole and Dorsal landmarks first!"
            return
        if self.spots is None:
            self.viewer.status = "Load spots first!"
            return

        self.viewer.status = "Orienting embryo..."
        pos = self.spots[["POSITION_X", "POSITION_Y", "POSITION_Z"]].values
        pos_new, params = orient_embryo(pos, self.animal_pole, self.dorsal_mark)

        self.spots["POSITION_X"] = pos_new[:, 0]
        self.spots["POSITION_Y"] = pos_new[:, 1]
        self.spots["POSITION_Z"] = pos_new[:, 2]
        self.transform_params = params
        self.oriented = True

        # Transform landmark positions into the new coordinate system
        mid = params["center"]
        R1, R2 = params["R1"], params["R2"]
        if self.animal_pole is not None:
            self.animal_pole = R2 @ (R1 @ (self.animal_pole - mid))
        if self.dorsal_mark is not None:
            self.dorsal_mark = R2 @ (R1 @ (self.dorsal_mark - mid))

        # Update layers
        self._refresh_points()
        self._refresh_tracks()
        self._update_landmarks()

        # Reset camera to centre on the reoriented data
        self.viewer.reset_view()
        self.viewer.status = (
            "Oriented: AP → top (-Y), Dorsal → right (+X). Camera reset."
        )

    # ── Sphere fitting ────────────────────────────────────────────────

    def _fit_sphere(self):
        """Fit sphere to the nuclei positions (background thread)."""
        from qtpy.QtCore import QTimer

        if self.spots is None:
            self.viewer.status = "Load spots first!"
            return

        x = self.spots["POSITION_X"].values.copy()
        y = self.spots["POSITION_Y"].values.copy()
        z = self.spots["POSITION_Z"].values.copy()

        self.lbl_sphere.value = "⏳ Fitting sphere..."
        self.viewer.status = "Fitting sphere in background..."

        self._pending_sphere = None
        self._pending_sphere_err = None

        def _bg():
            try:
                self._pending_sphere = fit_sphere(x, y, z)
            except Exception as e:
                self._pending_sphere_err = str(e)

        t = threading.Thread(target=_bg, daemon=True)
        t.start()

        def _poll():
            if t.is_alive():
                QTimer.singleShot(300, _poll)
                return
            if self._pending_sphere_err:
                self.lbl_sphere.value = f"Error: {self._pending_sphere_err}"
                self.viewer.status = "Sphere fit failed"
                return
            self._finish_sphere_fit(self._pending_sphere)

        QTimer.singleShot(300, _poll)

    def _finish_sphere_fit(self, result):
        """Finalise sphere fit on the main thread."""
        self.sphere_fit_result = result

        x = self.spots["POSITION_X"].values
        y = self.spots["POSITION_Y"].values
        z = self.spots["POSITION_Z"].values

        sph = compute_spherical_coords(
            np.column_stack([x, y, z]), result["center"], result["radius"]
        )
        for key, vals in sph.items():
            self.spots[key] = vals

        c = result["center"]
        self.lbl_sphere.value = (
            f"R={result['radius']:.0f}µm  RMSE={result['rmse']:.1f}µm  "
            f"Cov={result['coverage'] * 100:.0f}%\n"
            f"Center: ({c[0]:.0f}, {c[1]:.0f}, {c[2]:.0f})"
        )

        self._add_sphere_overlay()
        self.viewer.status = (
            f"Sphere fit: R={result['radius']:.0f}µm, RMSE={result['rmse']:.1f}µm"
        )

    def _add_sphere_overlay(self):
        """Add sphere cap overlay (only the portion near data)."""
        if self.sphere_fit_result is None:
            return

        c = self.sphere_fit_result["center"]
        R = self.sphere_fit_result["radius"]

        for name in ["Sphere mesh"]:
            try:
                self.viewer.layers.remove(name)
            except (ValueError, KeyError):
                pass

        pts = self.spots[["POSITION_X", "POSITION_Y", "POSITION_Z"]].values
        verts, faces = make_sphere_cap_mesh(c, R, pts)
        # napari Surface: vertices in (z, y, x) order for 3D
        verts_napari = verts[:, [2, 1, 0]]
        values = np.ones(len(verts))
        self.viewer.add_surface(
            (verts_napari, faces, values),
            name="Sphere mesh",
            opacity=0.25,
            colormap="cyan",
        )

    # ── Recolouring ───────────────────────────────────────────────────

    def _recolor(self, mode: str):
        """Recolour spots AND tracks by the chosen mode."""
        if self.spots_layer is None or self.spots is None:
            self.viewer.status = "Load spots first!"
            return

        self._current_color_mode = mode
        df = self.spots
        if mode == "frame":
            v = df["FRAME"].values.astype(float)
        elif mode == "depth" and "SPHERICAL_DEPTH" in df.columns:
            v = df["SPHERICAL_DEPTH"].values
        elif mode == "theta" and "THETA_DEG" in df.columns:
            v = df["THETA_DEG"].values
        elif mode == "phi" and "PHI_DEG" in df.columns:
            v = df["PHI_DEG"].values
        elif mode == "speed":
            v = self._compute_spot_speed()
        elif mode == "radial_vel" and "RADIAL_VELOCITY_SMOOTH" in df.columns:
            v = df["RADIAL_VELOCITY_SMOOTH"].values.copy()
            self._recolor_diverging(v, "Radial Velocity")
            self._rebuild_tracks()
            return
        elif mode == "ingressing" and "INGRESSING" in df.columns:
            self._recolor_binary("INGRESSING", "red", "Ingressing")
            self._rebuild_tracks()
            return
        elif mode == "roi" and "IN_ROI" in df.columns:
            self._recolor_binary("IN_ROI", "yellow", "ROI")
            self._rebuild_tracks()
            return
        else:
            self.viewer.status = (
                f"'{mode}' not available — fit sphere / apply ROI first."
            )
            return

        v_norm = (v - np.nanmin(v)) / max(np.nanmax(v) - np.nanmin(v), 1e-10)
        self.spots_layer.face_color = _viridis_colors(v_norm)
        self._rebuild_tracks()
        self.viewer.status = f"Coloured by {mode}"

    def _recolor_binary(self, column: str, highlight_color: str, label: str):
        """Recolour spots by a boolean column (IN_ROI, INGRESSING, etc)."""
        flags = self.spots[column].values.astype(bool)
        color_map = {
            "red": np.array([0.95, 0.1, 0.35, 1.0]),  # magenta-red for ingressing
            "yellow": np.array([1, 0.9, 0, 1]),
            "green": np.array([0, 0.8, 0, 1]),
        }
        hi = color_map.get(highlight_color, np.array([1, 0, 0, 1]))
        lo = np.array([0.5, 0.5, 0.5, 0.12])  # very dim grey for non-highlighted
        colors = np.where(flags[:, None], hi, lo)
        self.spots_layer.face_color = colors
        n = int(flags.sum())
        self.viewer.status = f"Coloured by {label}: {n:,} highlighted"

    def _recolor_diverging(self, values: np.ndarray, label: str):
        """Diverging colourmap: blue (negative) — white (0) — red (positive)."""
        v = values.copy()
        nan_mask = np.isnan(v)
        v[nan_mask] = 0.0
        vmax = max(np.abs(np.nanpercentile(values[~nan_mask], [2, 98])).max(), 1e-6)
        v_clip = np.clip(v / vmax, -1, 1)
        # Map [-1, 1] → RGBA: blue → white → red
        colors = np.zeros((len(v), 4))
        colors[:, 3] = 1.0
        neg = v_clip < 0
        pos = v_clip >= 0
        # Negative: blue (0,0,1) → white (1,1,1)
        t_neg = 1 + v_clip[neg]  # 0..1 (0=full blue, 1=white)
        colors[neg, 0] = t_neg
        colors[neg, 1] = t_neg
        colors[neg, 2] = 1.0
        # Positive: white (1,1,1) → red (1,0,0)
        t_pos = v_clip[pos]  # 0..1 (0=white, 1=full red)
        colors[pos, 0] = 1.0
        colors[pos, 1] = 1.0 - t_pos
        colors[pos, 2] = 1.0 - t_pos
        # NaN → grey
        colors[nan_mask] = [0.5, 0.5, 0.5, 0.3]
        self.spots_layer.face_color = colors
        self.viewer.status = f"Coloured by {label} (blue=inward, red=outward)"

    def _compute_spot_speed(self) -> np.ndarray:
        """Compute per-spot displacement speed (µm/frame) from track data."""
        df = self.spots
        speed = np.full(len(df), np.nan)
        if "TRACK_ID" not in df.columns:
            return speed
        # Sort by track and frame for diff
        sorted_idx = df.sort_values(["TRACK_ID", "FRAME"]).index
        tid = df.loc[sorted_idx, "TRACK_ID"].values
        x = df.loc[sorted_idx, "POSITION_X"].values
        y = df.loc[sorted_idx, "POSITION_Y"].values
        z = df.loc[sorted_idx, "POSITION_Z"].values
        f = df.loc[sorted_idx, "FRAME"].values

        # Compute displacement to previous point in same track
        dx = np.diff(x, prepend=np.nan)
        dy = np.diff(y, prepend=np.nan)
        dz = np.diff(z, prepend=np.nan)
        df_val = np.diff(f, prepend=-9999)
        same_track = np.diff(tid.astype(float), prepend=-9999) == 0
        dist = np.sqrt(dx**2 + dy**2 + dz**2)
        dt = np.maximum(df_val.astype(float), 1.0)
        sp = np.where(same_track & (df_val > 0), dist / dt, np.nan)
        speed[sorted_idx] = sp
        return speed

    # ── Ingression analysis ───────────────────────────────────────────

    def _compute_radial_velocity(self):
        """Compute per-spot radial velocity relative to the fitted sphere."""
        if self.sphere_fit_result is None:
            self.viewer.status = "Fit sphere first!"
            return
        if self.spots is None:
            return

        center = self.sphere_fit_result["center"]
        self.viewer.status = "Computing smoothed radial velocity..."
        self.spots = compute_radial_velocity(self.spots, center, smooth_window=5)

        rv = self.spots["RADIAL_VELOCITY_SMOOTH"].values
        valid = ~np.isnan(rv)
        if valid.any():
            med = np.nanmedian(rv)
            q10 = np.nanpercentile(rv, 10)
            q90 = np.nanpercentile(rv, 90)
            self.lbl_ingression.value = (
                f"V_r(smooth): med={med:.3f} µm/f  [{q10:.2f}, {q90:.2f}]"
            )
        else:
            self.lbl_ingression.value = "V_r: no track data"

        self.viewer.status = (
            "Radial velocity computed (smoothed, window=5). "
            "Use 'Auto-Detect Threshold' or 'Flag Ingressing Cells' next."
        )

    def _auto_detect_threshold(self):
        """Auto-detect the ingression V_r threshold using Otsu's method."""
        if self.spots is None or "TRACK_MEDIAN_RADIAL_VEL" not in self.spots.columns:
            self.viewer.status = "Compute radial velocity first!"
            return

        thresh, desc = auto_ingression_threshold(self.spots)
        # Clamp to slider range
        thresh = max(
            self.ingression_threshold.min, min(self.ingression_threshold.max, thresh)
        )
        self.ingression_threshold.value = thresh
        self.lbl_ingression.value = f"Auto: {desc}"
        self.viewer.status = f"Auto-threshold: {desc}"

    def _flag_ingression(self):
        """Flag ingressing cells based on multi-criteria detection."""
        if "RADIAL_VELOCITY" not in self.spots.columns:
            self.viewer.status = "Compute radial velocity first!"
            return

        thresh = self.ingression_threshold.value
        min_f = int(self.ingression_min_frames.value)
        inward_frac = self.ingression_inward_frac.value
        self.spots = flag_ingressing(
            self.spots,
            threshold=thresh,
            min_track_frames=min_f,
            min_inward_fraction=inward_frac,
        )

        n_ing = int(self.spots["INGRESSING"].sum())
        n_total = len(self.spots)
        n_tracks = 0
        if "TRACK_ID" in self.spots.columns:
            n_tracks = int(
                self.spots.loc[self.spots["INGRESSING"], "TRACK_ID"].nunique()
            )
        pct = n_ing / max(n_total, 1) * 100
        self.lbl_ingression.value = (
            f"Ingressing: {n_ing:,} spots ({n_tracks} tracks, {pct:.1f}%)\n"
            f"V_r≤{thresh:.3f} µm/f, inward≥{inward_frac:.0%}, ≥{min_f} frames"
        )
        self.viewer.status = (
            f"Flagged {n_ing:,} ingressing spots ({n_tracks} tracks). "
            f"Use 'Colour Ingressing' to visualise."
        )

    def _on_ingression_tracks_toggle(self):
        """Show only tracks/spots of ingressing cells, or all."""
        if self.spots is None:
            return
        if self.ingression_tracks_only.value and "INGRESSING" not in self.spots.columns:
            self.viewer.status = "Flag ingressing cells first!"
            self.ingression_tracks_only.value = False
            return
        self._refresh_points()  # also filter spots display
        self._rebuild_tracks()

    # ── Refresh layers ────────────────────────────────────────────────

    def _refresh_points(self):
        """Update the points layer data after orientation or slider change."""
        if self.spots_layer is None:
            return
        data, colors = self._build_display_points()
        self.spots_layer.data = data
        self.spots_layer.face_color = colors

    def _refresh_tracks(self):
        """Update the tracks layer after orientation or slider change."""
        self._rebuild_tracks()

    def _on_display_pct_changed(self):
        """Handle display percentage slider change."""
        self._display_pct = self.display_pct_slider.value
        if self.spots is not None:
            self._refresh_points()

    def _on_max_tracks_changed(self):
        """Handle max tracks slider change (debounced)."""
        self._max_tracks = int(self.max_tracks_slider.value)
        self._schedule_track_rebuild()

    # ── Track time range ──────────────────────────────────────────────

    def _on_track_frame_changed(self):
        """Handle track frame-range slider change (debounced)."""
        f_start = int(self.track_frame_start.value)
        f_end = int(self.track_frame_end.value)
        if f_end < f_start:
            f_end = f_start
        self._track_frame_min = f_start
        self._track_frame_max = f_end
        self._schedule_track_rebuild()

    def _schedule_track_rebuild(self):
        """Debounce track rebuilds — waits 200ms after last slider change."""
        from qtpy.QtCore import QTimer

        if not hasattr(self, "_debounce_timer"):
            self._debounce_timer = QTimer()
            self._debounce_timer.setSingleShot(True)
            self._debounce_timer.timeout.connect(self._deferred_track_rebuild)
        self._debounce_timer.start(200)  # ms

    def _deferred_track_rebuild(self):
        """Execute the actual track rebuild after debounce delay."""
        if self.spots is not None:
            self._rebuild_tracks()

    # ── 3D ROI box ────────────────────────────────────────────────────

    def _activate_roi_drawing(self):
        """Switch to 2D, add Shapes layer in rectangle-draw mode."""
        if self.spots is None:
            self.viewer.status = "Load spots first!"
            return

        self._drawing_roi = True  # suppress track rebuilds

        # Remove previous shapes layer
        if self._roi_shapes_layer is not None:
            try:
                self.viewer.layers.remove(self._roi_shapes_layer)
            except (ValueError, KeyError):
                pass
            self._roi_shapes_layer = None

        # Remember 3D state, switch to 2D (Shapes drawing only works in 2D)
        self._was_3d = self.viewer.dims.ndisplay == 3
        if self._was_3d:
            self.viewer.dims.ndisplay = 2

        self._roi_shapes_layer = self.viewer.add_shapes(
            name="ROI draw",
            ndim=4,
            edge_color="yellow",
            edge_width=2,
            face_color=np.array([1.0, 1.0, 0.0, 0.15]),
        )
        self._roi_shapes_layer.events.data.connect(self._on_roi_drawn)
        self._roi_shapes_layer.mode = "add_rectangle"
        self.viewer.status = (
            "\u25a3 Draw a rectangle to select a region "
            "(selects all nuclei across all Z depths; view returns to 3D after)"
        )

    def _on_roi_drawn(self, event=None):
        """Auto-apply ROI when a rectangle is drawn, then restore 3D."""
        if self._roi_updating:
            return
        if self._roi_shapes_layer is None:
            return
        shapes = self._roi_shapes_layer.data
        if len(shapes) == 0 or self.spots is None:
            return

        self._roi_updating = True
        try:
            # Keep only the last drawn rectangle
            if len(shapes) > 1:
                self._roi_shapes_layer.data = [shapes[-1]]
                return  # re-triggers via event; guard prevents recursion

            # In 2D mode with 4D data, rectangle verts are (4, 4): [t, z, y, x]
            rect_verts = np.asarray(shapes[0])
            if rect_verts.shape[1] == 4:
                x_vals = rect_verts[:, 3]  # POSITION_X
                y_vals = rect_verts[:, 2]  # POSITION_Y
            else:
                x_vals = rect_verts[:, 2]  # fallback 3D
                y_vals = rect_verts[:, 1]

            x_min, x_max = float(x_vals.min()), float(x_vals.max())
            y_min, y_max = float(y_vals.min()), float(y_vals.max())

            z_all = self.spots["POSITION_Z"].values
            z_min, z_max = float(z_all.min()), float(z_all.max())

            self._roi_bounds = {
                "x_min": x_min,
                "x_max": x_max,
                "y_min": y_min,
                "y_max": y_max,
                "z_min": z_min,
                "z_max": z_max,
            }

            in_roi = (
                (self.spots["POSITION_X"].values >= x_min)
                & (self.spots["POSITION_X"].values <= x_max)
                & (self.spots["POSITION_Y"].values >= y_min)
                & (self.spots["POSITION_Y"].values <= y_max)
            )
            self.spots["IN_ROI"] = in_roi

            n_roi = int(in_roi.sum())
            n_total = len(self.spots)
            self.lbl_roi.value = (
                f"ROI: {n_roi:,}/{n_total:,} ({n_roi / n_total * 100:.1f}%) selected"
            )

            # Remove the drawing layer and switch back to 3D
            try:
                self.viewer.layers.remove(self._roi_shapes_layer)
            except (ValueError, KeyError):
                pass
            self._roi_shapes_layer = None

            # Restore 3D view and draw wireframe
            if getattr(self, "_was_3d", True):
                self.viewer.dims.ndisplay = 3
            self._draw_roi_box()

            self.viewer.status = (
                f"ROI applied: {n_roi:,} spots inside rectangle (all Z depths)"
            )
        finally:
            self._roi_updating = False
            self._drawing_roi = False  # release guard — allow track rebuilds
            # Now do a single rebuild with the new ROI
            if self._roi_only:
                self._refresh_points()
                self._rebuild_tracks()

    def _draw_roi_box(self):
        """Draw the ROI bounding box as a wireframe Shapes layer."""
        if self._roi_bounds is None:
            return

        b = self._roi_bounds
        # Remove old
        for name in ["ROI box"]:
            try:
                self.viewer.layers.remove(name)
            except (ValueError, KeyError):
                pass

        # 12 edges of the box, in (z, y, x) for napari 3D
        corners = np.array(
            [
                [b["z_min"], b["y_min"], b["x_min"]],
                [b["z_min"], b["y_min"], b["x_max"]],
                [b["z_min"], b["y_max"], b["x_min"]],
                [b["z_min"], b["y_max"], b["x_max"]],
                [b["z_max"], b["y_min"], b["x_min"]],
                [b["z_max"], b["y_min"], b["x_max"]],
                [b["z_max"], b["y_max"], b["x_min"]],
                [b["z_max"], b["y_max"], b["x_max"]],
            ]
        )
        edges = [
            [0, 1],
            [0, 2],
            [1, 3],
            [2, 3],  # z_min face
            [4, 5],
            [4, 6],
            [5, 7],
            [6, 7],  # z_max face
            [0, 4],
            [1, 5],
            [2, 6],
            [3, 7],  # connecting edges
        ]
        lines = [np.array([corners[e[0]], corners[e[1]]]) for e in edges]

        self.viewer.add_shapes(
            lines,
            shape_type="line",
            edge_color="yellow",
            edge_width=3,
            name="ROI box",
            ndim=3,
        )

    def _on_roi_toggle(self):
        """Toggle: show all spots/tracks or only ROI-selected."""
        self._roi_only = self.roi_show_only.value
        if self.spots is not None:
            self._refresh_points()
            self._rebuild_tracks()

    def _clear_roi(self):
        """Clear the ROI selection, shapes layer, and wireframe box."""
        self._roi_bounds = None
        self._roi_only = False
        self._drawing_roi = False
        self.roi_show_only.value = False
        if self.spots is not None and "IN_ROI" in self.spots.columns:
            self.spots.drop(columns=["IN_ROI"], inplace=True)
        self.lbl_roi.value = "ROI: not set"
        for name in ["ROI box", "ROI draw"]:
            try:
                self.viewer.layers.remove(name)
            except (ValueError, KeyError):
                pass
        self._roi_shapes_layer = None
        if self.spots is not None:
            self._refresh_points()
            self._rebuild_tracks()
        self.viewer.status = "ROI cleared"

    # ── Export & Video ────────────────────────────────────────────────

    def _export(self):
        """Export comprehensive analysis output for downstream R pipeline."""
        import json

        out_dir = Path("analysis_output")
        out_dir.mkdir(exist_ok=True)
        files_written = []

        # 1. Enriched spots CSV (all computed columns)
        self.spots.to_csv(out_dir / "oriented_spots.csv", index=False)
        files_written.append("oriented_spots.csv")

        # 2. Sphere parameters
        if self.sphere_fit_result is not None:
            c = self.sphere_fit_result["center"]
            params = pd.DataFrame(
                {
                    "parameter": [
                        "center_x",
                        "center_y",
                        "center_z",
                        "radius",
                        "rmse",
                        "coverage",
                        "surface_used",
                        "cap_height",
                        "cap_base_radius",
                    ],
                    "value": [
                        c[0],
                        c[1],
                        c[2],
                        self.sphere_fit_result["radius"],
                        self.sphere_fit_result["rmse"],
                        self.sphere_fit_result["coverage"],
                        self.sphere_fit_result["surface_used"],
                        self.sphere_fit_result["cap_height"],
                        self.sphere_fit_result["cap_base_radius"],
                    ],
                }
            )
            params.to_csv(out_dir / "sphere_params.csv", index=False)
            files_written.append("sphere_params.csv")

        # 3. ROI bounds
        if self._roi_bounds is not None:
            roi_df = pd.DataFrame(
                {
                    "parameter": list(self._roi_bounds.keys()),
                    "value": list(self._roi_bounds.values()),
                }
            )
            roi_df.to_csv(out_dir / "roi_bounds.csv", index=False)
            files_written.append("roi_bounds.csv")

        # 4. Ingression parameters
        if "INGRESSING" in self.spots.columns:
            n_ing = int(self.spots["INGRESSING"].sum())
            n_tracks_ing = 0
            if "TRACK_ID" in self.spots.columns:
                n_tracks_ing = int(
                    self.spots.loc[
                        self.spots["INGRESSING"].astype(bool), "TRACK_ID"
                    ].nunique()
                )
            ing_df = pd.DataFrame(
                {
                    "parameter": [
                        "threshold_um_per_frame",
                        "min_track_frames",
                        "min_inward_fraction",
                        "detection_method",
                        "n_ingressing_spots",
                        "n_ingressing_tracks",
                    ],
                    "value": [
                        self.ingression_threshold.value,
                        int(self.ingression_min_frames.value),
                        self.ingression_inward_frac.value,
                        "multi-criteria (smoothed median V_r + inward fraction)",
                        n_ing,
                        n_tracks_ing,
                    ],
                }
            )
            ing_df.to_csv(out_dir / "ingression_params.csv", index=False)
            files_written.append("ingression_params.csv")

        # 5. Analysis metadata (JSON)
        meta = {
            "oriented": self.oriented,
            "sphere_fitted": self.sphere_fit_result is not None,
            "roi_set": self._roi_bounds is not None,
            "time_window": {
                "frame_min": self._track_frame_min,
                "frame_max": self._track_frame_max,
            },
            "display": {
                "display_pct": self._display_pct,
                "max_tracks": self._max_tracks,
                "current_color_mode": self._current_color_mode,
            },
            "n_spots": len(self.spots),
            "n_frames": int(self.spots["FRAME"].nunique()),
            "columns_exported": list(self.spots.columns),
        }
        if self.sphere_fit_result is not None:
            meta["sphere_radius"] = float(self.sphere_fit_result["radius"])
        if self._roi_bounds is not None:
            meta["roi_bounds"] = {k: float(v) for k, v in self._roi_bounds.items()}
        with open(out_dir / "analysis_metadata.json", "w") as f:
            json.dump(meta, f, indent=2, default=str)
        files_written.append("analysis_metadata.json")

        self.lbl_export.value = f"Exported {len(files_written)} files to {out_dir}/"
        self.viewer.status = f"Exported: {', '.join(files_written)}"

    def _record_video(self):
        """Record a time-lapse video iterating through the time window."""
        from qtpy.QtCore import QTimer

        if self.spots is None:
            self.viewer.status = "Load spots first!"
            return

        out_dir = Path("analysis_output")
        out_dir.mkdir(exist_ok=True)
        f_start = int(self.track_frame_start.value)
        f_end = int(self.track_frame_end.value)
        all_frames = sorted(self.spots["FRAME"].unique())
        frames = [int(f) for f in all_frames if f_start <= f <= f_end]
        if not frames:
            self.viewer.status = "No frames in selected range!"
            return

        self._vid_frames_list = frames
        self._vid_images = []
        self._vid_idx = 0
        self._vid_total = len(frames)
        self._vid_path = out_dir / "embryo_timelapse.mp4"
        self.lbl_export.value = f"Recording: 0/{self._vid_total}..."
        self.viewer.status = "Recording video \u2014 please wait..."

        def _capture_next():
            if self._vid_idx >= self._vid_total:
                # --- write video ---
                fps = max(1, min(30, self._vid_total // 10))
                try:
                    import imageio

                    writer = imageio.get_writer(
                        str(self._vid_path),
                        fps=fps,
                        codec="libx264",
                        quality=8,
                        macro_block_size=1,
                    )
                    for img in self._vid_images:
                        writer.append_data(img)
                    writer.close()
                except Exception:
                    try:
                        import imageio

                        writer = imageio.get_writer(
                            str(self._vid_path),
                            fps=fps,
                            macro_block_size=1,
                        )
                        for img in self._vid_images:
                            writer.append_data(img)
                        writer.close()
                    except Exception as e2:
                        # Fallback: save PNGs
                        png_dir = out_dir / "video_frames"
                        png_dir.mkdir(exist_ok=True)
                        import imageio

                        for i, img in enumerate(self._vid_images):
                            imageio.imwrite(str(png_dir / f"frame_{i:04d}.png"), img)
                        self.lbl_export.value = f"Saved {len(self._vid_images)} PNGs"
                        self.viewer.status = (
                            f"mp4 failed ({e2}), saved PNGs to {png_dir}/"
                        )
                        return
                self.lbl_export.value = f"Video: {self._vid_path}"
                self.viewer.status = (
                    f"Video saved to {self._vid_path} ({self._vid_total} frames)"
                )
                return

            frame = self._vid_frames_list[self._vid_idx]
            dims = list(self.viewer.dims.current_step)
            dims[0] = frame
            self.viewer.dims.current_step = tuple(dims)

            def _screenshot():
                img = self.viewer.screenshot(canvas_only=True)
                self._vid_images.append(img)
                self._vid_idx += 1
                if self._vid_idx % 10 == 0:
                    self.lbl_export.value = (
                        f"Recording: {self._vid_idx}/{self._vid_total}..."
                    )
                QTimer.singleShot(50, _capture_next)

            QTimer.singleShot(150, _screenshot)

        QTimer.singleShot(200, _capture_next)


# ############################################################################
#                               UTILITIES
# ############################################################################


def _viridis_colors(values: np.ndarray) -> np.ndarray:
    """Map normalised [0,1] values to RGBA viridis colours."""
    try:
        import matplotlib.pyplot as plt

        cmap = plt.cm.viridis
        return cmap(values)
    except ImportError:
        pass

    # Fallback: simple blue→green→yellow gradient
    n = len(values)
    colors = np.zeros((n, 4))
    colors[:, 0] = np.clip(values * 2 - 1, 0, 1)  # R
    colors[:, 1] = np.clip(1 - 2 * np.abs(values - 0.5), 0, 1)  # G
    colors[:, 2] = np.clip(1 - values * 2, 0, 1)  # B
    colors[:, 3] = 1.0  # A
    return colors


# ############################################################################
#                                  MAIN
# ############################################################################


def main():
    parser = argparse.ArgumentParser(
        description="napari-based 4D embryo viewer for TrackMate data"
    )
    parser.add_argument(
        "spots",
        nargs="?",
        default=None,
        help="Path to spots CSV (optional — can load from UI)",
    )
    parser.add_argument(
        "--tracks", "-t", default=None, help="Path to tracks CSV (optional)"
    )
    parser.add_argument(
        "--max-spots",
        type=int,
        default=None,
        help="Max spots to load (subsample for speed)",
    )
    args = parser.parse_args()

    spots_df = None
    tracks_df = None

    if args.spots:
        print(f"Loading spots from {args.spots}...")
        spots_df = load_spots(args.spots)
        print(f"  {len(spots_df):,} spots, {spots_df['FRAME'].nunique()} frames")

        if args.max_spots and len(spots_df) > args.max_spots:
            print(f"  Subsampling to {args.max_spots:,} spots...")
            spots_df = spots_df.sample(args.max_spots, random_state=42)

        if args.tracks:
            print(f"Loading tracks from {args.tracks}...")
            tracks_df = load_trackmate_csv(args.tracks)
            print(f"  {len(tracks_df):,} track records")
    else:
        print("Launching viewer — use 'Load Spots CSV' button to open files.")

    print("Launching napari viewer...")
    viewer = EmbryoViewer(spots_df, tracks_df)

    napari, _ = _import_napari()
    napari.run()


if __name__ == "__main__":
    main()
