# =============================================================================
# napari-based 4D Embryo Viewer — Orientation, Sphere Fit & Ingression Analysis
# =============================================================================
#
# PURPOSE:
#   Interactive 4D (3D + time) visualisation and annotation of light-sheet
#   microscopy nuclear-tracking data (TrackMate or ultrack exports) for medaka
#   and zebrafish gastrulation analysis.
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
    os.environ["QT_QPA_PLATFORM"] = "wayland"
os.environ.setdefault("QT_AUTO_SCREEN_SCALE_FACTOR", "1")
os.environ.setdefault("QT_FONT_DPI", "96")
# Suppress harmless Qt C++ warnings (QSocketNotifier, DPI messages)
os.environ.setdefault("QT_LOGGING_RULES", "*.warning=false")

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


def load_ultrack_csv(path: str | Path) -> pd.DataFrame:
    """Load an ultrack tracks CSV and map columns to internal format.

    Ultrack exports columns: track_id, t, z, y, x, id, parent_track_id, parent_id.
    Maps to: TRACK_ID, FRAME, POSITION_Z, POSITION_Y, POSITION_X.
    """
    path = Path(path)
    df = pd.read_csv(path, low_memory=False)

    col_map = {
        "track_id": "TRACK_ID",
        "t": "FRAME",
        "z": "POSITION_Z",
        "y": "POSITION_Y",
        "x": "POSITION_X",
    }
    missing = [c for c in col_map if c not in df.columns]
    if missing:
        raise ValueError(f"Not a valid ultrack CSV — missing columns: {missing}")

    df = df.rename(columns=col_map)
    for c in ["POSITION_X", "POSITION_Y", "POSITION_Z", "FRAME", "TRACK_ID"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["FRAME"] = df["FRAME"].astype(int)

    # Keep parent info if present
    if "parent_track_id" in df.columns:
        df["PARENT_TRACK_ID"] = pd.to_numeric(df["parent_track_id"], errors="coerce")

    return df


def _detect_csv_format(path: str | Path) -> str:
    """Detect whether a CSV is ultrack or TrackMate format.

    Returns 'ultrack' or 'trackmate'.
    """
    path = Path(path)
    df_head = pd.read_csv(path, nrows=5, low_memory=False)
    ultrack_cols = {"track_id", "t", "z", "y", "x"}
    if ultrack_cols.issubset(set(df_head.columns)):
        return "ultrack"
    return "trackmate"


def load_spots(path: str | Path) -> pd.DataFrame:
    """Load spots CSV and return cleaned DataFrame.

    Auto-detects ultrack vs TrackMate format.
    """
    fmt = _detect_csv_format(path)
    if fmt == "ultrack":
        df = load_ultrack_csv(path)
    else:
        df = load_trackmate_csv(path)

    required = ["POSITION_X", "POSITION_Y", "POSITION_Z", "FRAME"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    df = df.dropna(subset=required)
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
    # Use np.maximum to avoid divide-by-zero; results where dt<=0 are
    # masked to NaN by the np.where condition anyway.
    dt_safe = np.maximum(dt, 1e-12)
    vel = np.where(same_track & (dt > 0), dr / dt_safe, np.nan)

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
    """Flag tracks as ingressing using embryology-informed multi-criteria detection.

    During gastrulation in medaka/zebrafish, ingressing cells move
    radially inward through the blastoderm, transitioning from the
    epiblast to the hypoblast (or mesendoderm). True ingression is
    characterised by:
      - Sustained inward radial motion (not transient fluctuation)
      - Net radial displacement toward the sphere centre
      - Consistently negative radial velocity over most of the track

    Criteria — a track is flagged INGRESSING if ALL of:
      1. Track has ≥ *min_track_frames* time points
      2. Track **median** smoothed V_r ≤ *threshold*
      3. ≥ *min_inward_fraction* of track points have V_r < 0
      4. Net radial displacement is inward (first→last distance decreases)
      5. Track has ≥ 3 consecutive frames of inward movement (sustained)

    Also computes per-track columns for downstream analysis:
      INGRESSION_ONSET_FRAME  – first frame of sustained inward movement
      TRACK_NET_RADIAL_DISP   – net change in radial distance (neg = inward)
      TRACK_SUSTAINED_INWARD  – max consecutive inward frames

    Parameters
    ----------
    df : spots DataFrame with TRACK_MEDIAN_RADIAL_VEL, TRACK_INWARD_FRACTION
    threshold : µm/frame – tracks with median V_r ≤ threshold → ingressing
    min_track_frames : minimum track length
    min_inward_fraction : minimum fraction of inward-moving points

    Returns
    -------
    DataFrame with INGRESSING (bool) column + enrichment columns
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

    # ── Net radial displacement criterion ──
    # True ingressing cells end up closer to the sphere centre than
    # they started (epiblast→hypoblast transition).
    net_disp = np.full(len(df), np.nan)
    sustained_inward = np.zeros(len(df), dtype=int)
    onset_frame = np.full(len(df), np.nan)

    if "RADIAL_DIST_TO_CENTER" in df.columns:
        sorted_df = df.sort_values(["TRACK_ID", "FRAME"])

        # Vectorised net radial displacement: last - first per track
        first_r = sorted_df.groupby("TRACK_ID")["RADIAL_DIST_TO_CENTER"].transform(
            "first"
        )
        last_r = sorted_df.groupby("TRACK_ID")["RADIAL_DIST_TO_CENTER"].transform(
            "last"
        )
        net_per_spot = (last_r - first_r).values
        net_disp[sorted_df.index] = net_per_spot

        # Sustained inward: max consecutive frames with V_r < 0 (per track)
        # This requires run-length logic so we iterate per track
        if "RADIAL_VELOCITY_SMOOTH" in sorted_df.columns:
            for tid, grp in sorted_df.groupby("TRACK_ID"):
                idx = grp.index
                vr = grp["RADIAL_VELOCITY_SMOOTH"].values
                max_run, cur_run, first_inward = 0, 0, np.nan
                for j, v in enumerate(vr):
                    if not np.isnan(v) and v < 0:
                        cur_run += 1
                        if cur_run == 1 and np.isnan(first_inward):
                            first_inward = grp["FRAME"].values[j]
                        if cur_run > max_run:
                            max_run = cur_run
                    else:
                        if cur_run < 3:  # reset onset if run was too short
                            first_inward = np.nan
                        cur_run = 0
                sustained_inward[idx] = max_run
                onset_frame[idx] = first_inward

    df["TRACK_NET_RADIAL_DISP"] = net_disp
    df["TRACK_SUSTAINED_INWARD"] = sustained_inward
    df["INGRESSION_ONSET_FRAME"] = onset_frame

    # Net displacement is inward — require meaningful movement, not noise
    net_inward = np.where(
        np.isnan(net_disp),
        True,  # if no data, don't block
        net_disp < -0.5,  # at least 0.5 µm net inward
    )

    # Sustained movement: scale with track length — at least 15% of frames
    # or 3, whichever is larger.  Short noise bursts get filtered out while
    # long tracks need proportionally more consecutive inward evidence.
    min_sustained = np.maximum(3, (track_len * 0.15).astype(int))
    sustained_ok = sustained_inward >= min_sustained

    df["INGRESSING"] = (
        long_enough
        & (~np.isnan(vel))
        & (vel <= threshold)
        & inward_ok
        & net_inward
        & sustained_ok
    )

    # ── Continuous ingression score (0–1) for visualisation ──
    # Combines all criteria into a smooth score so the visualisation can
    # show confidence rather than binary on/off.
    abs_thresh = max(abs(threshold), 0.01)
    s_vel = np.clip((threshold - vel) / abs_thresh, 0, 2) / 2
    s_vel = np.nan_to_num(s_vel, nan=0.0)

    if "TRACK_INWARD_FRACTION" in df.columns:
        s_inward = np.clip(df["TRACK_INWARD_FRACTION"].values, 0, 1)
    else:
        s_inward = np.full(len(df), 0.5)

    s_sustained = np.clip(sustained_inward / np.maximum(track_len.values, 1), 0, 1)

    disp_for_score = np.where(np.isnan(net_disp), 0.0, -net_disp)
    s_disp = np.clip(
        disp_for_score / np.maximum(track_len.values.astype(float), 1), 0, 1
    )

    raw_score = 0.30 * s_vel + 0.20 * s_inward + 0.25 * s_sustained + 0.25 * s_disp
    # Zero out score for non-ingressing tracks to keep it interpretable
    raw_score[~df["INGRESSING"].values] = 0.0
    df["INGRESSION_SCORE"] = raw_score

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
#                       CROSS-SECTION VIEWER
# ############################################################################


class CrossSectionViewer:
    """A separate window showing a 2D cross-section through the 3D point cloud.

    The user selects a cutting axis (X, Y, or Z) and a position along that
    axis.  All spots within a configurable slab thickness around that position
    are projected onto the remaining two axes as a 2D scatter, colored by
    the selected colour mode.

    Mirrors the main viewer's controls: colour modes, display %, tracks
    within the slab, sphere intersection, time slider.
    """

    # Mapping from cutting axis to the 2 display axes
    _AXES_MAP = {
        "X": ("POSITION_Z", "POSITION_Y"),  # (vertical, horizontal)
        "Y": ("POSITION_Z", "POSITION_X"),
        "Z": ("POSITION_Y", "POSITION_X"),
    }
    _AXES_LABELS = {
        "X": ("Z", "Y"),
        "Y": ("Z", "X"),
        "Z": ("Y", "X"),
    }
    # Maps POSITION_* column names to xyz index for coordinate lookups
    _COL_TO_IDX = {"POSITION_X": 0, "POSITION_Y": 1, "POSITION_Z": 2}

    def __init__(self, parent: "EmbryoViewer"):
        napari, _ = _import_napari()
        from magicgui.widgets import (
            Container,
            FloatSlider,
            Label,
            ComboBox,
            PushButton,
            CheckBox,
        )

        self.parent = parent
        self.df = parent.spots
        self._axis = "Y"  # cutting axis
        self._thickness = 20.0  # slab total thickness
        self._color_mode = parent._current_color_mode
        self._display_pct = 100.0
        self._show_all = False
        self._track_width = 2.0
        self._show_tracks = True

        # Set max tracks to actual track count
        if "TRACK_ID" in self.df.columns:
            n_total_tracks = int(self.df["TRACK_ID"].nunique())
        else:
            n_total_tracks = 500
        self._max_tracks = n_total_tracks

        # Compute axis ranges
        self._ranges = {
            "X": (
                float(self.df["POSITION_X"].min()),
                float(self.df["POSITION_X"].max()),
            ),
            "Y": (
                float(self.df["POSITION_Y"].min()),
                float(self.df["POSITION_Y"].max()),
            ),
            "Z": (
                float(self.df["POSITION_Z"].min()),
                float(self.df["POSITION_Z"].max()),
            ),
        }
        fmin = int(self.df["FRAME"].min())
        fmax = int(self.df["FRAME"].max())

        # Separate napari viewer in 2D
        self.viewer = napari.Viewer(title="Cross-Section View", ndisplay=2)
        self.points_layer = None
        self.tracks_layer = None
        self._sphere_layer = None
        self._landmarks_layer = None
        self._center_layer = None
        self._recording = False
        self._saved_view_state = None  # camera + layer settings snapshot

        r = self._ranges[self._axis]
        mid = (r[0] + r[1]) / 2.0

        # ── Slab controls ──
        self.axis_combo = ComboBox(value="Y", choices=["X", "Y", "Z"], label="Cut axis")
        self.axis_combo.changed.connect(self._on_axis_changed)

        self.pos_slider = FloatSlider(
            value=mid, min=r[0], max=r[1], step=1.0, label="Position"
        )
        self.pos_slider.changed.connect(self._schedule_rebuild)

        self.thickness_slider = FloatSlider(
            value=20.0, min=1.0, max=200.0, step=1.0, label="Slab thickness"
        )
        self.thickness_slider.changed.connect(self._on_thickness_changed)

        # Quick-jump buttons for common positions
        btn_jump_center = PushButton(text="Jump to sphere center")
        btn_jump_center.changed.connect(self._jump_to_center)
        btn_jump_ap = PushButton(text="Jump to animal pole")
        btn_jump_ap.changed.connect(self._jump_to_ap)

        # ── Time ──
        self.frame_slider = FloatSlider(
            value=fmin, min=fmin, max=fmax, step=1, label="Frame"
        )
        self.frame_slider.changed.connect(self._schedule_rebuild)

        self.all_frames_check = CheckBox(value=False, text="Show all frames")
        self.all_frames_check.changed.connect(self._on_all_frames_changed)

        # ── Display ──
        self.display_pct_slider = FloatSlider(
            value=100, min=1, max=100, step=1, label="Display %"
        )
        self.display_pct_slider.changed.connect(self._on_display_pct_changed)

        self.show_tracks_check = CheckBox(value=True, text="Show tracks")
        self.show_tracks_check.changed.connect(self._on_show_tracks_changed)

        self.max_tracks_slider = FloatSlider(
            value=n_total_tracks,
            min=100,
            max=max(500, n_total_tracks),
            step=500,
            label="Max tracks",
        )
        self.max_tracks_slider.changed.connect(self._on_max_tracks_changed)

        self.track_width_slider = FloatSlider(
            value=2.0, min=0.5, max=10.0, step=0.5, label="Track width"
        )
        self.track_width_slider.changed.connect(self._on_track_width_changed)

        self.point_size_slider = FloatSlider(
            value=3.0, min=1.0, max=20.0, step=0.5, label="Point size"
        )
        self.point_size_slider.changed.connect(self._schedule_rebuild)

        # ── Colour mode buttons ──
        btn_color_frame = PushButton(text="Colour by Frame")
        btn_color_frame.changed.connect(lambda: self._set_color("frame"))
        btn_color_depth = PushButton(text="Colour by Depth")
        btn_color_depth.changed.connect(lambda: self._set_color("depth"))
        btn_color_depth_travel = PushButton(text="Colour by Depth Travel")
        btn_color_depth_travel.changed.connect(lambda: self._set_color("depth_travel"))
        btn_color_theta = PushButton(text="Colour by Latitude")
        btn_color_theta.changed.connect(lambda: self._set_color("theta"))
        btn_color_phi = PushButton(text="Colour by Longitude")
        btn_color_phi.changed.connect(lambda: self._set_color("phi"))
        btn_color_speed = PushButton(text="Colour by Speed")
        btn_color_speed.changed.connect(lambda: self._set_color("speed"))
        btn_color_radvel = PushButton(text="Colour by Radial Velocity")
        btn_color_radvel.changed.connect(lambda: self._set_color("radial_vel"))
        btn_color_track = PushButton(text="Colour by Track ID")
        btn_color_track.changed.connect(lambda: self._set_color("track_id"))

        # ── Annotation toggles ──
        self.show_landmarks_check = CheckBox(
            value=True, text="Show landmarks (AP/Dorsal)"
        )
        self.show_landmarks_check.changed.connect(self._schedule_rebuild)
        self.show_sphere_check = CheckBox(value=True, text="Show sphere intersection")
        self.show_sphere_check.changed.connect(self._schedule_rebuild)
        self.show_center_check = CheckBox(value=True, text="Show sphere center")
        self.show_center_check.changed.connect(self._schedule_rebuild)

        # ── Record video ──
        self.frame_start_slider = FloatSlider(
            value=fmin, min=fmin, max=fmax, step=1, label="Start frame"
        )
        self.frame_end_slider = FloatSlider(
            value=fmax, min=fmin, max=fmax, step=1, label="End frame"
        )
        self.frame_step_slider = FloatSlider(
            value=1, min=1, max=10, step=1, label="Frame step"
        )
        btn_record = PushButton(text="\U0001f3ac Record Time-lapse Video")
        btn_record.changed.connect(self._record_video)
        btn_screenshot = PushButton(text="\U0001f4f7 Screenshot")
        btn_screenshot.changed.connect(self._save_screenshot)

        self.lbl_info = Label(value="")
        self.lbl_annotations = Label(value="")
        self.lbl_record = Label(value="")

        container = Container(
            widgets=[
                Label(value="── Slab ──"),
                self.axis_combo,
                self.pos_slider,
                self.thickness_slider,
                btn_jump_center,
                btn_jump_ap,
                Label(value="── Time ──"),
                self.frame_slider,
                self.all_frames_check,
                Label(value="── Display ──"),
                self.display_pct_slider,
                self.point_size_slider,
                self.show_tracks_check,
                self.max_tracks_slider,
                self.track_width_slider,
                Label(value="── Annotations ──"),
                self.show_landmarks_check,
                self.show_sphere_check,
                self.show_center_check,
                Label(value="── Colour ──"),
                btn_color_frame,
                btn_color_depth,
                btn_color_depth_travel,
                btn_color_theta,
                btn_color_phi,
                btn_color_speed,
                btn_color_radvel,
                btn_color_track,
                Label(value="── Record ──"),
                self.frame_start_slider,
                self.frame_end_slider,
                self.frame_step_slider,
                btn_record,
                btn_screenshot,
                self.lbl_record,
                Label(value="── Info ──"),
                self.lbl_info,
                self.lbl_annotations,
            ]
        )

        from qtpy.QtWidgets import QScrollArea
        from qtpy.QtCore import Qt

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setMinimumWidth(280)
        scroll.setWidget(container.native)
        self.viewer.window.add_dock_widget(scroll, name="Controls", area="right")

        self._rebuild()

    # ── Quick-jump helpers ───────────────────────────────────────────

    def _jump_to_center(self):
        """Move the slice position to the sphere center along the current axis."""
        sph = self.parent.sphere_fit_result
        if sph is None:
            self.viewer.status = "Fit sphere in main viewer first!"
            return
        axis_idx = {"X": 0, "Y": 1, "Z": 2}[self._axis]
        val = sph["center"][axis_idx]
        self.pos_slider.value = np.clip(val, self.pos_slider.min, self.pos_slider.max)
        self.viewer.status = f"Jumped to sphere center along {self._axis}"

    def _jump_to_ap(self):
        """Move the slice position to the animal pole along the current axis."""
        ap = self.parent.animal_pole
        if ap is None:
            self.viewer.status = "Set animal pole in main viewer first!"
            return
        axis_idx = {"X": 0, "Y": 1, "Z": 2}[self._axis]
        val = ap[axis_idx]
        self.pos_slider.value = np.clip(val, self.pos_slider.min, self.pos_slider.max)
        self.viewer.status = f"Jumped to animal pole along {self._axis}"

    # ── Coordinate helpers ───────────────────────────────────────────

    def _project_point_xyz(self, xyz: np.ndarray) -> np.ndarray:
        """Project an (x, y, z) world point to 2D [v, h] for current axis."""
        v_col, h_col = self._AXES_MAP[self._axis]
        v = xyz[self._COL_TO_IDX[v_col]]
        h = xyz[self._COL_TO_IDX[h_col]]
        return np.array([v, h])

    def _axis_value_of_point(self, xyz: np.ndarray) -> float:
        """Return the coordinate of xyz along the current cutting axis."""
        axis_idx = {"X": 0, "Y": 1, "Z": 2}[self._axis]
        return float(xyz[axis_idx])

    # ── Widget callbacks ─────────────────────────────────────────────

    def _on_axis_changed(self):
        self._axis = self.axis_combo.value
        r = self._ranges[self._axis]
        mid = (r[0] + r[1]) / 2.0
        self.pos_slider.min = r[0]
        self.pos_slider.max = r[1]
        self.pos_slider.value = mid
        self._rebuild()

    def _on_thickness_changed(self):
        self._thickness = self.thickness_slider.value
        self._schedule_rebuild()

    def _on_all_frames_changed(self):
        self._show_all = self.all_frames_check.value
        self._rebuild()

    def _on_display_pct_changed(self):
        self._display_pct = self.display_pct_slider.value
        self._rebuild()

    def _on_show_tracks_changed(self):
        self._show_tracks = self.show_tracks_check.value
        self._rebuild()

    def _on_max_tracks_changed(self):
        self._max_tracks = int(self.max_tracks_slider.value)
        self._schedule_rebuild()

    def _on_track_width_changed(self):
        self._track_width = float(self.track_width_slider.value)
        self._schedule_rebuild()

    def _set_color(self, mode: str):
        needed = {
            "depth": "SPHERICAL_DEPTH",
            "depth_travel": "SPHERICAL_DEPTH",
            "theta": "THETA_DEG",
            "phi": "PHI_DEG",
            "radial_vel": "RADIAL_VELOCITY_SMOOTH",
        }
        if mode in needed and needed[mode] not in self.df.columns:
            self.viewer.status = f"'{mode}' not available — fit sphere first."
            return
        self._color_mode = mode
        self._rebuild()

    def _schedule_rebuild(self):
        """Debounce rapid slider changes."""
        from qtpy.QtCore import QTimer

        if not hasattr(self, "_debounce_timer"):
            self._debounce_timer = QTimer()
            self._debounce_timer.setSingleShot(True)
            self._debounce_timer.timeout.connect(self._rebuild)
        self._debounce_timer.start(150)

    # ── Core rebuild ─────────────────────────────────────────────────

    def _get_slab_mask(self) -> np.ndarray:
        """Return boolean mask for spots within the current slab + frame."""
        df = self.df
        axis_col = f"POSITION_{self._axis}"
        pos = self.pos_slider.value
        half = self._thickness / 2.0
        mask = (df[axis_col].values >= pos - half) & (df[axis_col].values <= pos + half)

        if not self._show_all:
            frame = int(self.frame_slider.value)
            frame_mask = df["FRAME"].values == frame
            mask = mask & frame_mask
        return mask

    def _project(self, idx: np.ndarray) -> np.ndarray:
        """Project indexed spots to 2D [v, h] coordinates."""
        v_col, h_col = self._AXES_MAP[self._axis]
        h = self.df[h_col].values[idx]
        v = self.df[v_col].values[idx]
        return np.column_stack([v, h])

    def _rebuild(self):
        """Full rebuild of points, tracks, annotations, and sphere layers."""
        # Remove all managed layers
        for attr in (
            "points_layer",
            "tracks_layer",
            "_sphere_layer",
            "_landmarks_layer",
            "_center_layer",
        ):
            layer = getattr(self, attr, None)
            if layer is not None:
                try:
                    self.viewer.layers.remove(layer)
                except (ValueError, KeyError):
                    pass
                setattr(self, attr, None)

        mask = self._get_slab_mask()
        idx = np.where(mask)[0]

        # Display % subsampling
        if self._display_pct < 100 and len(idx) > 0:
            n_show = max(100, int(len(idx) * self._display_pct / 100))
            if n_show < len(idx):
                rng = np.random.default_rng(0)
                idx = rng.choice(idx, n_show, replace=False)
                idx.sort()

        # ── Points ──
        if len(idx) > 0:
            data = self._project(idx)
            colors = self._get_colors(idx)
            self.points_layer = self.viewer.add_points(
                data,
                name="Nuclei",
                size=self.point_size_slider.value,
                face_color=colors,
                border_width=0,
                ndim=2,
            )

        # ── Tracks within slab ──
        if self._show_tracks and len(idx) > 0:
            self._add_slab_tracks(mask)

        # ── Sphere cross-section circle ──
        if self.show_sphere_check.value:
            self._add_sphere_intersection()

        # ── Sphere center marker ──
        if self.show_center_check.value:
            self._add_center_marker()

        # ── Landmark markers (AP / Dorsal) ──
        if self.show_landmarks_check.value:
            self._add_landmark_markers()

        # ── Axis labels ──
        v_label, h_label = self._AXES_LABELS[self._axis]
        self.viewer.dims.axis_labels = [v_label, h_label]

        # ── Info label ──
        pos = self.pos_slider.value
        half = self._thickness / 2.0
        info = f"Cut {self._axis}={pos:.0f} ± {half:.0f}  |  {len(idx):,} spots"
        if not self._show_all:
            info += f"  |  frame {int(self.frame_slider.value)}"
        else:
            info += "  |  all frames"
        n_tracks_shown = 0
        if self.tracks_layer is not None:
            try:
                n_tracks_shown = len(np.unique(self.tracks_layer.data[:, 0]))
            except Exception:
                pass
        if n_tracks_shown > 0:
            info += f"  |  {n_tracks_shown:,} tracks"
        self.lbl_info.value = info

        # ── Annotations info ──
        ann_parts = []
        sph = self.parent.sphere_fit_result
        if sph is not None:
            c = sph["center"]
            ann_parts.append(
                f"Sphere: R={sph['radius']:.0f}, center=({c[0]:.0f},{c[1]:.0f},{c[2]:.0f})"
            )
        if self.parent.animal_pole is not None:
            ap = self.parent.animal_pole
            ann_parts.append(f"AP: ({ap[0]:.0f},{ap[1]:.0f},{ap[2]:.0f})")
        if self.parent.dorsal_mark is not None:
            dm = self.parent.dorsal_mark
            ann_parts.append(f"Dorsal: ({dm[0]:.0f},{dm[1]:.0f},{dm[2]:.0f})")
        self.lbl_annotations.value = (
            "\\n".join(ann_parts) if ann_parts else "No annotations set"
        )

        if not self._recording:
            self.viewer.reset_view()
        else:
            self._restore_view_state()

    def _snapshot_view_state(self):
        """Save camera and all layer display settings before recording."""
        cam = self.viewer.camera
        self._saved_view_state = {
            "camera": {
                "center": tuple(cam.center),
                "zoom": cam.zoom,
                "angles": tuple(cam.angles),
            },
            "layers": {},
        }
        for layer in self.viewer.layers:
            props = {
                "visible": layer.visible,
                "opacity": layer.opacity,
            }
            # Tracks-specific settings
            if hasattr(layer, "tail_length"):
                props["tail_length"] = layer.tail_length
            if hasattr(layer, "tail_width"):
                props["tail_width"] = layer.tail_width
            # Points-specific settings
            if hasattr(layer, "size") and not callable(layer.size):
                try:
                    props["size"] = (
                        float(layer.size)
                        if np.ndim(layer.size) == 0
                        else layer.size.copy()
                    )
                except Exception:
                    pass
            if hasattr(layer, "symbol"):
                props["symbol"] = layer.symbol
            # Shapes-specific settings
            if hasattr(layer, "edge_width"):
                props["edge_width"] = layer.edge_width
            self._saved_view_state["layers"][layer.name] = props

    def _restore_view_state(self):
        """Restore camera and layer display settings from snapshot."""
        if self._saved_view_state is None:
            return
        # Restore camera
        cam_state = self._saved_view_state["camera"]
        self.viewer.camera.center = cam_state["center"]
        self.viewer.camera.zoom = cam_state["zoom"]
        # Restore layer settings
        saved_layers = self._saved_view_state["layers"]
        for layer in self.viewer.layers:
            if layer.name in saved_layers:
                props = saved_layers[layer.name]
                for key, val in props.items():
                    try:
                        setattr(layer, key, val)
                    except Exception:
                        pass

    def _add_slab_tracks(self, slab_mask: np.ndarray):
        """Add 2D tracks for track segments within the slab."""
        df = self.df
        if "TRACK_ID" not in df.columns:
            return

        # Get track IDs that have spots in slab
        slab_tids = df.loc[slab_mask, "TRACK_ID"].dropna().unique()
        if len(slab_tids) == 0:
            return

        # Subsample if needed
        if len(slab_tids) > self._max_tracks:
            rng = np.random.default_rng(42)
            slab_tids = rng.choice(slab_tids, self._max_tracks, replace=False)

        # Get ALL points of those tracks (not just slab points) so
        # trajectories are continuous, then filter by slab + time
        tid_set = set(slab_tids)
        in_tracks = df["TRACK_ID"].isin(tid_set).values
        if not self._show_all:
            # For single-frame mode, show full track history of those tracks
            # within the slab for context
            axis_col = f"POSITION_{self._axis}"
            pos = self.pos_slider.value
            half = self._thickness / 2.0
            in_slab = (df[axis_col].values >= pos - half) & (
                df[axis_col].values <= pos + half
            )
            track_mask = in_tracks & in_slab
        else:
            track_mask = in_tracks & slab_mask

        df_t = df.loc[track_mask].copy()
        if df_t.empty:
            return

        track_len = df_t.groupby("TRACK_ID")["FRAME"].transform("count")
        df_t = df_t[track_len >= 2]
        if df_t.empty:
            return

        df_t = df_t.sort_values(["TRACK_ID", "FRAME"])

        v_col, h_col = self._AXES_MAP[self._axis]
        tracks_arr = np.column_stack(
            [
                df_t["TRACK_ID"].values,
                df_t["FRAME"].values,
                df_t[v_col].values,
                df_t[h_col].values,
            ]
        )

        # Color tracks consistently with points
        props = {}
        color_by = "track_id"
        cmap = "hsv"

        col_map = {
            "depth": "SPHERICAL_DEPTH",
            "theta": "THETA_DEG",
            "phi": "PHI_DEG",
            "radial_vel": "RADIAL_VELOCITY_SMOOTH",
            "frame": "FRAME",
        }
        mode = self._color_mode
        if mode == "depth_travel" and "SPHERICAL_DEPTH" in df_t.columns:
            depth_range = df_t.groupby("TRACK_ID")["SPHERICAL_DEPTH"].transform(
                lambda s: s.max() - s.min()
            ).values.astype(float)
            nan_m = np.isnan(depth_range)
            if not nan_m.all():
                depth_range[nan_m] = 0.0
                valid = depth_range[~nan_m]
                if len(valid) > 0:
                    vmax = max(np.percentile(valid, 98), 1e-6)
                    depth_range = np.clip(depth_range / vmax, 0, 1)
                props["color_value"] = depth_range
                color_by = "color_value"
                cmap = "plasma"
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
            tail_width=self._track_width,
            tail_length=40,
            color_by=color_by,
            colormap=cmap,
        )
        if props:
            kw["properties"] = props

        try:
            self.tracks_layer = self.viewer.add_tracks(tracks_arr, **kw)
        except Exception:
            pass  # tracks can fail if too few points; not critical

    def _add_sphere_intersection(self):
        """Draw the sphere's intersection circle with the cutting plane."""
        sph = self.parent.sphere_fit_result
        if sph is None:
            return

        center = sph["center"]
        R = sph["radius"]
        pos = self.pos_slider.value

        axis_idx = {"X": 0, "Y": 1, "Z": 2}[self._axis]
        d = abs(pos - center[axis_idx])
        if d >= R:
            return  # plane doesn't intersect sphere

        r = np.sqrt(R**2 - d**2)
        cv, ch = self._project_point_xyz(center)

        theta = np.linspace(0, 2 * np.pi, 200)
        circle_v = cv + r * np.sin(theta)
        circle_h = ch + r * np.cos(theta)
        circle_pts = np.column_stack([circle_v, circle_h])

        self._sphere_layer = self.viewer.add_shapes(
            [circle_pts],
            shape_type="polygon",
            edge_color="cyan",
            edge_width=2,
            face_color=np.array([0.0, 1.0, 1.0, 0.05]),
            name="Sphere section",
            ndim=2,
        )

    def _add_center_marker(self):
        """Draw a crosshair at the sphere center projected into the slice."""
        sph = self.parent.sphere_fit_result
        if sph is None:
            return
        center_2d = self._project_point_xyz(sph["center"])
        arm = sph["radius"] * 0.05  # 5% of radius
        lines = [
            np.array(
                [[center_2d[0] - arm, center_2d[1]], [center_2d[0] + arm, center_2d[1]]]
            ),
            np.array(
                [[center_2d[0], center_2d[1] - arm], [center_2d[0], center_2d[1] + arm]]
            ),
        ]
        self._center_layer = self.viewer.add_shapes(
            lines,
            shape_type="line",
            edge_color="white",
            edge_width=2,
            name="Sphere center",
            ndim=2,
        )

    def _add_landmark_markers(self):
        """Show animal pole and dorsal landmarks projected into the slice."""
        pts = []
        colors = []
        texts = []

        for name, xyz, color in [
            ("AP", self.parent.animal_pole, "red"),
            ("Dorsal", self.parent.dorsal_mark, "cyan"),
        ]:
            if xyz is None:
                continue
            projected = self._project_point_xyz(xyz)
            pts.append(projected)
            texts.append(name)
            colors.append(color)

        if not pts:
            return

        pts_arr = np.array(pts)
        self._landmarks_layer = self.viewer.add_points(
            pts_arr,
            name="Landmarks",
            size=25,
            face_color=colors,
            symbol="diamond",
            border_width=1,
            border_color="white",
            ndim=2,
            text={
                "string": texts,
                "color": "white",
                "anchor": "upper_left",
                "size": 14,
                "translation": [5, 5],
            },
        )

    def _get_colors(self, idx: np.ndarray) -> np.ndarray:
        """Color the cross-section points."""
        if len(idx) == 0:
            return np.array([[0.5, 0.5, 0.5, 1.0]])

        df = self.df
        mode = self._color_mode

        # Track ID coloring: random per-track colors
        if mode == "track_id" and "TRACK_ID" in df.columns:
            tids = df["TRACK_ID"].values[idx]
            unique_tids = np.unique(tids[~np.isnan(tids)])
            if len(unique_tids) > 0:
                try:
                    import matplotlib.pyplot as plt

                    cmap = plt.cm.tab20
                except ImportError:
                    cmap = None
                tid_to_color = {}
                for i, tid in enumerate(unique_tids):
                    if cmap is not None:
                        tid_to_color[tid] = cmap(i % 20)
                    else:
                        rng = np.random.default_rng(int(tid))
                        tid_to_color[tid] = (*rng.random(3), 1.0)
                colors = np.array(
                    [tid_to_color.get(t, (0.5, 0.5, 0.5, 0.3)) for t in tids]
                )
                return colors

        # Diverging mode for radial velocity
        if mode == "radial_vel" and "RADIAL_VELOCITY_SMOOTH" in df.columns:
            v = df["RADIAL_VELOCITY_SMOOTH"].values[idx].copy()
            nan_mask = np.isnan(v)
            v[nan_mask] = 0.0
            valid = v[~nan_mask]
            vmax = (
                max(np.abs(np.nanpercentile(valid, [2, 98])).max(), 1e-6)
                if len(valid) > 0
                else 1.0
            )
            v_clip = np.clip(v / vmax, -1, 1)
            colors = np.zeros((len(v), 4))
            colors[:, 3] = 1.0
            neg = v_clip < 0
            pos = v_clip >= 0
            t_neg = 1 + v_clip[neg]
            colors[neg, 0] = t_neg
            colors[neg, 1] = t_neg
            colors[neg, 2] = 1.0
            t_pos = v_clip[pos]
            colors[pos, 0] = 1.0
            colors[pos, 1] = 1.0 - t_pos
            colors[pos, 2] = 1.0 - t_pos
            colors[nan_mask] = [0.5, 0.5, 0.5, 0.3]
            return colors

        # Speed mode
        if mode == "speed":
            speed = self._compute_speed(idx)
            v = speed.astype(float)
        elif mode == "depth" and "SPHERICAL_DEPTH" in df.columns:
            v = df["SPHERICAL_DEPTH"].values[idx]
        elif mode == "depth_travel" and "SPHERICAL_DEPTH" in df.columns and "TRACK_ID" in df.columns:
            depth_range = df.groupby("TRACK_ID")["SPHERICAL_DEPTH"].transform(
                lambda s: s.max() - s.min()
            )
            v = depth_range.values[idx].astype(float)
            nan_mask = np.isnan(v)
            v[nan_mask] = 0.0
            valid = v[~nan_mask]
            if len(valid) > 0:
                vmax = max(np.percentile(valid, 98), 1e-6)
                v_norm = np.clip(v / vmax, 0, 1)
            else:
                v_norm = np.zeros_like(v)
            colors = _viridis_colors(v_norm)
            colors[nan_mask] = [0.4, 0.4, 0.4, 0.15]
            return colors
        elif mode == "theta" and "THETA_DEG" in df.columns:
            v = df["THETA_DEG"].values[idx]
        elif mode == "phi" and "PHI_DEG" in df.columns:
            v = df["PHI_DEG"].values[idx]
        else:
            v = df["FRAME"].values[idx].astype(float)

        v = v.astype(float)
        vmin, vmax = np.nanmin(v), np.nanmax(v)
        v_norm = (v - vmin) / max(vmax - vmin, 1e-10)
        return _viridis_colors(v_norm)

    def _save_screenshot(self):
        """Save a single screenshot of the current cross-section view."""
        out_dir = Path("analysis_output")
        out_dir.mkdir(exist_ok=True)
        axis = self._axis
        pos = int(self.pos_slider.value)
        frame = int(self.frame_slider.value)
        fname = f"cross_section_{axis}_{pos}_f{frame}.png"
        path = out_dir / fname
        img = self.viewer.screenshot(canvas_only=True)
        try:
            import imageio

            imageio.imwrite(str(path), img)
            self.lbl_record.value = f"Saved: {path}"
            self.viewer.status = f"Screenshot saved to {path}"
        except Exception as e:
            self.lbl_record.value = f"Error: {e}"

    def _record_video(self):
        """Record a time-lapse video — same approach as main embryo viewer.

        Builds all slab data with a time dimension (ndim=3), then steps
        through frames with dims.current_step and takes screenshots.
        No manual per-frame layer manipulation needed.
        """
        from qtpy.QtCore import QTimer

        f_start = int(self.frame_start_slider.value)
        f_end = int(self.frame_end_slider.value)
        f_step = max(1, int(self.frame_step_slider.value))

        all_frames = sorted(self.df["FRAME"].unique())
        frames = [int(f) for f in all_frames if f_start <= f <= f_end][::f_step]
        if not frames:
            self.viewer.status = "No frames in selected range!"
            return

        # ── Snapshot current state ──
        self._snapshot_view_state()
        self._recording = True

        # ── Remove all existing managed layers ──
        for attr in (
            "points_layer",
            "tracks_layer",
            "_sphere_layer",
            "_landmarks_layer",
            "_center_layer",
        ):
            layer = getattr(self, attr, None)
            if layer is not None:
                try:
                    self.viewer.layers.remove(layer)
                except (ValueError, KeyError):
                    pass
                setattr(self, attr, None)

        # ── Build 3D (time, v, h) points for ALL slab frames ──
        df = self.df
        axis_col = f"POSITION_{self._axis}"
        pos = self.pos_slider.value
        half = self._thickness / 2.0
        slab_mask = (df[axis_col].values >= pos - half) & (
            df[axis_col].values <= pos + half
        )

        v_col, h_col = self._AXES_MAP[self._axis]
        slab_idx = np.where(slab_mask)[0]

        # Display % subsampling per frame
        if self._display_pct < 100 and len(slab_idx) > 0:
            keep = np.zeros(len(df), dtype=bool)
            for f in frames:
                fm = slab_mask & (df["FRAME"].values == f)
                fidx = np.where(fm)[0]
                n_show = max(100, int(len(fidx) * self._display_pct / 100))
                if n_show < len(fidx):
                    rng = np.random.default_rng(0)
                    fidx = rng.choice(fidx, n_show, replace=False)
                keep[fidx] = True
            slab_idx = np.where(slab_mask & keep)[0]

        saved = self._saved_view_state or {}
        saved_layers = saved.get("layers", {})
        nuclei_visible = saved_layers.get("Nuclei", {}).get("visible", True)
        tracks_visible = saved_layers.get("Tracks", {}).get(
            "visible", self._show_tracks
        )

        if nuclei_visible and len(slab_idx) > 0:
            # 3D data: [frame, v, h]
            data_3d = np.column_stack(
                [
                    df["FRAME"].values[slab_idx],
                    df[v_col].values[slab_idx],
                    df[h_col].values[slab_idx],
                ]
            )
            colors = self._get_colors(slab_idx)
            pt_size = saved_layers.get("Nuclei", {}).get(
                "size", self.point_size_slider.value
            )
            try:
                pt_size = float(pt_size)
            except Exception:
                pt_size = self.point_size_slider.value
            self.points_layer = self.viewer.add_points(
                data_3d,
                name="Nuclei",
                size=pt_size,
                face_color=colors,
                border_width=0,
                ndim=3,
            )

        # ── Build tracks (already has time column — works with 3D) ──
        if tracks_visible and self._show_tracks and "TRACK_ID" in df.columns:
            self._add_slab_tracks(slab_mask)
            if self.tracks_layer is not None:
                tprops = saved_layers.get("Tracks", {})
                for key in ("tail_length", "tail_width", "opacity"):
                    if key in tprops:
                        try:
                            setattr(self.tracks_layer, key, tprops[key])
                        except Exception:
                            pass

        # ── Re-add annotations as 3D (visible at all times) ──
        if self.show_sphere_check.value:
            self._add_sphere_intersection()
        if self.show_center_check.value:
            self._add_center_marker()
        if self.show_landmarks_check.value:
            self._add_landmark_markers()

        # ── Restore camera (adding 3D layers may have changed it) ──
        self._restore_view_state()

        out_dir = Path("analysis_output")
        out_dir.mkdir(exist_ok=True)
        axis = self._axis
        slab_pos = int(pos)
        vid_path = out_dir / f"cross_section_{axis}_{slab_pos}.mp4"

        self._vid_frames = frames
        self._vid_images = []
        self._vid_idx = 0
        self._vid_total = len(frames)
        self._vid_path = vid_path
        self.lbl_record.value = f"Recording: 0/{self._vid_total}..."
        self.viewer.status = "Recording cross-section video \u2014 please wait..."

        def _capture_next():
            if self._vid_idx >= self._vid_total:
                self._finish_video()
                return

            frame = self._vid_frames[self._vid_idx]
            # Step through time — exactly like the main embryo viewer
            dims = list(self.viewer.dims.current_step)
            dims[0] = frame
            self.viewer.dims.current_step = tuple(dims)

            def _screenshot():
                try:
                    img = self.viewer.screenshot(canvas_only=True)
                except (RuntimeError, AttributeError):
                    # Widget was closed during recording
                    return
                if img.ndim == 3 and img.shape[2] == 4:
                    img = img[:, :, :3]
                h, w = img.shape[:2]
                if h % 2 or w % 2:
                    img = img[: h - h % 2, : w - w % 2]
                self._vid_images.append(img)
                self._vid_idx += 1
                if self._vid_idx % 10 == 0 or self._vid_idx == self._vid_total:
                    self.lbl_record.value = (
                        f"Recording: {self._vid_idx}/{self._vid_total}..."
                    )
                QTimer.singleShot(30, _capture_next)

            QTimer.singleShot(150, _screenshot)

        QTimer.singleShot(200, _capture_next)

    def _finish_video(self):
        """Write captured frames to mp4 (or fallback to PNGs), then restore 2D view."""
        fps = max(1, min(30, self._vid_total // 10))
        out_dir = Path("analysis_output")
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
            self.lbl_record.value = f"Video: {self._vid_path}"
            self.viewer.status = (
                f"Video saved to {self._vid_path} ({self._vid_total} frames, {fps} fps)"
            )
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
                self.lbl_record.value = f"Video: {self._vid_path}"
                self.viewer.status = f"Video saved to {self._vid_path}"
            except Exception as e2:
                png_dir = out_dir / "cross_section_frames"
                png_dir.mkdir(exist_ok=True)
                import imageio

                for i, img in enumerate(self._vid_images):
                    imageio.imwrite(str(png_dir / f"frame_{i:04d}.png"), img)
                self.lbl_record.value = f"Saved {len(self._vid_images)} PNGs"
                self.viewer.status = f"mp4 failed ({e2}), saved PNGs to {png_dir}/"

        # Restore 2D view: set frame slider to last recorded frame,
        # exit recording mode, do a normal 2D rebuild.
        if self._vid_frames:
            last_frame = self._vid_frames[-1]
            self.frame_slider.value = last_frame
        self._recording = False
        self._saved_view_state = None
        self._rebuild()

    def _compute_speed(self, idx: np.ndarray) -> np.ndarray:
        """Per-spot displacement speed for selected indices."""
        df = self.df
        speed = np.full(len(df), np.nan)
        if "TRACK_ID" not in df.columns:
            return speed[idx]
        sorted_i = df.sort_values(["TRACK_ID", "FRAME"]).index
        tid = df.loc[sorted_i, "TRACK_ID"].values
        x = df.loc[sorted_i, "POSITION_X"].values
        y = df.loc[sorted_i, "POSITION_Y"].values
        z = df.loc[sorted_i, "POSITION_Z"].values
        f = df.loc[sorted_i, "FRAME"].values.astype(float)
        dx = np.diff(x, prepend=np.nan)
        dy = np.diff(y, prepend=np.nan)
        dz = np.diff(z, prepend=np.nan)
        df_v = np.diff(f, prepend=-9999)
        same = np.diff(tid.astype(float), prepend=-9999) == 0
        dist = np.sqrt(dx**2 + dy**2 + dz**2)
        dt = np.maximum(df_v, 1.0)
        sp = np.where(same & (df_v > 0), dist / dt, np.nan)
        speed[sorted_i] = sp
        return speed[idx]


# ############################################################################
#                            NAPARI VIEWER
# ############################################################################


class EmbryoViewer:
    """napari-based 4D embryo viewer with orientation & sphere fitting."""

    def __init__(
        self,
        spots_df: pd.DataFrame | None = None,
        tracks_df: pd.DataFrame | None = None,
        segments_zarr_path: str | Path | None = None,
        voxel_size: tuple[float, float, float] = (1.0, 1.0, 1.0),
        downsample: int = 1,
        preload: bool = False,
    ):
        napari, magicgui = _import_napari()

        self.spots = spots_df
        self.tracks_raw = tracks_df

        # State
        self.animal_pole: np.ndarray | None = None
        self.dorsal_mark: np.ndarray | None = None
        self.margin_landmark: np.ndarray | None = None       # margin Y cut
        self.ingression_landmark: np.ndarray | None = None   # ingression centre
        self.oriented = False
        self.sphere_fit_result: dict | None = None
        self.transform_params: dict | None = None
        self.spots_layer = None
        self.tracks_layer = None
        self._tracks_layer_secondary = None  # for dual-layer ingression view
        self._display_pct: float = 100.0
        self._max_tracks: int = 500000
        self._track_frame_min: int | None = None
        self._track_frame_max: int | None = None
        self._roi_bounds: dict | None = None  # {x_min, x_max, y_min, ...}
        self._roi_only: bool = False  # show only ROI spots
        self._roi_shapes_layer = None
        self._roi_updating: bool = False  # guard for event recursion
        self._drawing_roi: bool = False  # guard: suppress rebuilds while drawing
        self._current_color_mode: str = "frame"  # tracks colour mode
        self._tracks_cache_key: tuple | None = None  # for caching track builds
        self._display_idx: np.ndarray | None = (
            None  # indices into self.spots for displayed subset
        )
        self._track_width: float = 2.0  # track line width

        # Segments state
        self._segments_data = None  # dask array (T, Z, Y, X)
        self._segments_layer = None  # napari Labels layer
        self._voxel_size = voxel_size  # (Z, Y, X) in µm
        self._downsample_factor = downsample
        self._preload = preload

        # Create viewer
        self.viewer = napari.Viewer(title="Embryo 4D Viewer", ndisplay=3)

        # Load segments if provided
        if segments_zarr_path is not None:
            self._load_segments(segments_zarr_path)

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

    # ── Segments (Labels) layer ────────────────────────────────────────

    def _load_segments(self, zarr_path: str | Path) -> None:
        """Load segments.zarr and add a Labels layer.

        Two modes controlled by self._preload:

        **Lazy (default)** — keeps data on disk via dask.  Each timepoint
        is read on demand when the slider moves.  Works on any machine
        (needs only ~200 MB per displayed frame).  Scrubbing is slower
        because of disk I/O.

        **Preload (--preload)** — loads everything into RAM, remaps
        int32→uint16, and builds a GPU-optimized display volume for
        instant, buttery-smooth time scrubbing.  Keeps full-res data
        in RAM for analysis and a spatially-downsampled display copy
        (controlled by --downsample, default 2) for fast rendering.
        """
        import dask.array as da
        import zarr
        import time as _time

        zarr_path = Path(zarr_path).resolve()
        print(f"Loading segments: {zarr_path}")
        store = zarr.open(str(zarr_path), mode="r")

        if isinstance(store, zarr.Array):
            data = da.from_zarr(str(zarr_path))
            print(f"  4D zarr array: shape={data.shape}, dtype={data.dtype}")
        elif isinstance(store, zarr.Group):
            keys = sorted(int(k) for k in store.keys())
            first = store[str(keys[0])]
            lazy_frames = [
                da.from_zarr(str(zarr_path), component=str(k)) for k in keys
            ]
            data = da.stack(lazy_frames, axis=0)
            print(f"  Zarr group: {len(keys)} timepoints, "
                  f"spatial={first.shape}, dtype={first.dtype}")
        else:
            raise ValueError(f"Unexpected zarr layout at {zarr_path}")

        ds = self._downsample_factor

        if self._preload:
            # ── PRELOAD MODE: fast load, consistent colors, instant scrub ─
            import multiprocessing as mp
            from concurrent.futures import ThreadPoolExecutor
            n_cores = mp.cpu_count() or 1
            nbytes_raw = data.nbytes
            print(f"  Preloading {nbytes_raw / 1e9:.1f} GB into RAM "
                  f"({n_cores} cores) ...")
            t0 = _time.time()
            data_np = data.compute(scheduler="threads",
                                   num_workers=n_cores)
            print(f"  Loaded in {_time.time() - t0:.0f}s")
            del data

            # Keep full-resolution data for analysis
            self._segments_data = data_np

            # ── Build display frames ─────────────────────────────────
            # Each label gets a CONSISTENT colour across all timepoints
            # via Knuth multiplicative hash → uint8 intensity.
            # Same nucleus = same label ID = same colour always.
            #
            # Then: 3D Image layer + vispy direct texture swap bypasses
            # napari's entire 4D pipeline → zero Python overhead per
            # frame change.
            print(f"  Building {data_np.shape[0]} display frames "
                  f"({ds}× downsample) ...", end=" ", flush=True)
            t0 = _time.time()

            if ds > 1:
                display = data_np[:, ::ds, ::ds, ::ds]
            else:
                display = data_np

            # Hash-based global colour mapping (label → uint8 [5..255])
            # Knuth multiplicative hash gives good spread; background 0→0.
            def _hash_frame(t):
                vol = display[t]
                # Vectorised hash: (label * 2654435761) mod 251 + 5
                # Background (0) stays 0 → transparent.
                out = np.zeros(vol.shape, dtype=np.uint8)
                mask = vol != 0
                ids = vol[mask].astype(np.uint64)
                out[mask] = ((ids * np.uint64(2654435761)) % np.uint64(251)
                             + np.uint64(5)).astype(np.uint8)
                return np.ascontiguousarray(out)

            n_tp = display.shape[0]
            with ThreadPoolExecutor(max_workers=n_cores) as pool:
                self._segment_frames = list(pool.map(_hash_frame, range(n_tp)))

            frame_mb = self._segment_frames[0].nbytes / 1e6
            total_gb = n_tp * frame_mb / 1e3
            print(f"{n_tp} × {frame_mb:.1f} MB = "
                  f"{total_gb:.1f} GB ({_time.time() - t0:.1f}s)")

            # Transparent-background turbo colormap
            from napari.utils.colormaps import Colormap
            import matplotlib.pyplot as plt
            _lut = plt.get_cmap("turbo")(np.linspace(0, 1, 256))
            _lut[0] = [0, 0, 0, 0]  # background → fully transparent
            cmap = Colormap(colors=_lut.astype(np.float32),
                            name="turbo_nuclei")

            # Add as 3D layer — avoids napari's 4D slicing pipeline
            scale_3d = (ds, ds, ds)
            print(f"  Adding 3D Image layer ...")
            self._segments_layer = self.viewer.add_image(
                self._segment_frames[0],
                name="Nuclei 3D",
                scale=scale_3d,
                opacity=0.9,
                rendering="attenuated_mip",
                colormap=cmap,
                contrast_limits=[0, 255],
                blending="translucent",
            )

            # Kill thumbnail (scipy.ndimage.zoom on every refresh)
            self._segments_layer._update_thumbnail = lambda: None

            # ── CRITICAL PERFORMANCE FIX ─────────────────────────────
            # napari's viewer._update_layers is connected to
            # dims.events.current_step AND dims.events.point.
            # On every frame tick it calls _slice_dims on EVERY layer:
            #   → Points: np.min/max on 2M×4 array (extent) + slice filter
            #   → Tracks: full _set_view_slice + thumbnail + extent
            #   → Image:  slice + thumbnail + extent
            # This is the ENTIRE source of lag.  GPU sits idle waiting.
            #
            # Fix: disconnect napari's _update_layers from the time
            # slider events.  Our own minimal callback swaps only the
            # segments vispy texture.  Points/Tracks become static
            # overlays showing all timepoints at once (which is perfect
            # for trajectory visualisation).
            viewer_model = self.viewer
            try:
                viewer_model.dims.events.current_step.disconnect(
                    viewer_model._update_layers
                )
                viewer_model.dims.events.point.disconnect(
                    viewer_model._update_layers
                )
                print("  ✓ Disconnected napari per-frame layer processing "
                      "(eliminates 2M-point filtering per frame)")
            except (ValueError, TypeError) as exc:
                print(f"  ⚠ Could not disconnect _update_layers: {exc}")

            # Find vispy Volume node for direct GPU texture swap
            self._vispy_volume = None
            self._vispy_canvas = None
            try:
                qt_viewer = self.viewer.window._qt_viewer
                canvas = qt_viewer.canvas
                if hasattr(canvas, 'layer_to_visual'):
                    vl = canvas.layer_to_visual.get(self._segments_layer)
                    if vl is not None and hasattr(vl, 'node'):
                        node = vl.node
                        if hasattr(node, 'set_data'):
                            self._vispy_volume = node
                            self._vispy_canvas = canvas
                            print("  ✓ Vispy direct texture path active")
            except Exception as exc:
                print(f"  ⚠ Vispy path unavailable ({exc})")

            # Connect time slider → instant frame swap (ONLY segments)
            def _on_segments_time(event,
                                  frames=self._segment_frames,
                                  layer=self._segments_layer,
                                  vol=self._vispy_volume,
                                  canv=self._vispy_canvas):
                t = int(self.viewer.dims.current_step[0])
                if 0 <= t < len(frames):
                    if vol is not None:
                        vol.set_data(frames[t])
                        canv.native.update()
                    else:
                        layer.data = frames[t]

            self._segments_time_cb = _on_segments_time
            self.viewer.dims.events.current_step.connect(_on_segments_time)

            print(f"  Ready: {n_tp} tp  |  "
                  f"full-res {data_np.shape[1:]} "
                  f"({data_np.nbytes / 1e9:.1f} GB)  |  "
                  f"display {self._segment_frames[0].shape} "
                  f"({total_gb:.1f} GB)  |  "
                  f"scale {scale_3d}")
        else:
            # ── LAZY MODE: dask, low RAM, slower scrubbing ──
            if ds > 1:
                data = data[:, ::ds, ::ds, ::ds]
                print(f"  Lazy downsample {ds}× → effective shape {data.shape}")

            self._segments_data = data
            # Match tracks coordinate system (raw voxel units)
            scale_4d = (1, ds, ds, ds)

            self._segments_layer = self.viewer.add_labels(
                data,
                name="Nuclei 3D",
                scale=scale_4d,
                opacity=0.5,
            )

            est_frame_mb = (data[0].nbytes / 1e6) if hasattr(data[0], 'nbytes') else 0
            print(f"  Labels layer (lazy): {data.shape[0]} tp, "
                  f"~{est_frame_mb:.0f} MB/frame, scale={scale_4d}")

        self.viewer.dims.axis_labels = ["t", "z", "y", "x"]

    @staticmethod
    def _remap_4d_labels(data: np.ndarray) -> np.ndarray:
        """Remap a 4D int32/int64 label array to contiguous uint16 per frame.

        Parallelised across CPU cores using shared-memory output array.
        Avoids np.unique (which sorts the entire volume) — instead uses
        a scatter approach with np.flatnonzero on a boolean mask, which
        is 5-10× faster for sparse label volumes.
        """
        import multiprocessing as mp
        from concurrent.futures import ThreadPoolExecutor

        n_tp = data.shape[0]
        out = np.zeros(data.shape, dtype=np.uint16)

        def _remap_one(t: int) -> None:
            vol = data[t]
            if vol.max() == 0:
                return
            # Fast unique via boolean scatter — avoids sorting the
            # entire volume.  For ~5000 labels in 58M voxels this is
            # much faster than np.unique which sorts all 58M values.
            max_id = int(vol.max())
            seen = np.zeros(max_id + 1, dtype=bool)
            seen[vol.ravel()] = True
            seen[0] = False  # background
            unique_ids = np.flatnonzero(seen)
            n_labels = len(unique_ids)
            if n_labels == 0:
                return
            if n_labels > 65535:
                unique_ids = unique_ids[:65535]
                n_labels = 65535
            lut = np.zeros(max_id + 1, dtype=np.uint16)
            lut[unique_ids] = np.arange(1, n_labels + 1, dtype=np.uint16)
            out[t] = lut[vol]

        n_workers = min(n_tp, mp.cpu_count() or 1)
        # Use threads — the work is numpy C-level so releases the GIL
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            list(pool.map(_remap_one, range(n_tp)))

        return out

    def _load_segments_from_ui(self) -> None:
        """Open a folder dialog to load a segments.zarr directory."""
        from qtpy.QtWidgets import QFileDialog

        zarr_dir = QFileDialog.getExistingDirectory(
            self.viewer.window._qt_window,
            "Select segments.zarr directory",
        )
        if not zarr_dir:
            return
        try:
            # Remove existing segments layer + callback
            if self._segments_layer is not None:
                # Disconnect time callback if exists
                cb = getattr(self, '_segments_time_cb', None)
                if cb is not None:
                    try:
                        self.viewer.dims.events.current_step.disconnect(cb)
                    except (ValueError, TypeError):
                        pass
                    self._segments_time_cb = None
                try:
                    self.viewer.layers.remove(self._segments_layer)
                except (ValueError, KeyError):
                    pass
                self._segments_layer = None
            self._load_segments(zarr_dir)
            self.lbl_segments_status.value = (
                f"Loaded: {self._segments_data.shape}"
            )
        except Exception as e:
            self.lbl_segments_status.value = f"Error: {e}"

    def _on_segments_opacity_changed(self, value: float) -> None:
        if self._segments_layer is not None:
            self._segments_layer.opacity = value / 100.0

    def _on_segments_visible_changed(self, value: bool) -> None:
        if self._segments_layer is not None:
            self._segments_layer.visible = value

    def _add_spots_layer(self):
        """Add the main nuclei points layer."""
        data, idx = self._build_display_points()
        self._display_idx = idx
        colors = self._get_display_colors()
        self.spots_layer = self.viewer.add_points(
            data,
            name="Nuclei",
            size=2,
            face_color=colors,
            border_width=0,
            ndim=4,
            out_of_slice_display=False,
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
        return data, idx

    def _get_display_colors(self) -> np.ndarray:
        """Compute face colors for the currently displayed subset based on current color mode."""
        idx = self._display_idx
        if idx is None or self.spots is None:
            return np.array([[0.5, 0.5, 0.5, 1.0]])
        df = self.spots
        mode = self._current_color_mode

        # ── Ingression mode: color by inward radial velocity ──
        # Each spot is colored by how strongly it is moving inward at
        # that moment.  Non-ingressing spots are ghostly grey; ingressing
        # spots transition from warm orange (gentle inward) to bright
        # magenta (strong inward dive), making the ingression movement
        # itself visually pop.
        if mode == "ingressing" and "INGRESSING" in df.columns:
            flags = df["INGRESSING"].values[idx].astype(bool)
            n_pts = len(idx)
            colors = np.tile([0.45, 0.45, 0.45, 0.06], (n_pts, 1))
            if flags.any() and "RADIAL_VELOCITY_SMOOTH" in df.columns:
                vr = df["RADIAL_VELOCITY_SMOOTH"].values[idx][flags].copy()
                vr = np.nan_to_num(vr, nan=0.0)
                # Normalise: 0 = no inward motion, 1 = strongest inward
                # Only negative V_r matters; clip positives to 0
                inward = np.clip(-vr, 0, None)  # flip sign: positive = inward
                vmax = (
                    np.percentile(inward[inward > 0], 95) if (inward > 0).any() else 1.0
                )
                t = np.clip(inward / max(vmax, 1e-6), 0, 1)
                # Gradient: warm orange (t=0) → hot magenta (t=0.5) → bright white-pink (t=1)
                colors[flags, 0] = 0.85 + 0.15 * t  # R: stays high
                colors[flags, 1] = 0.35 * (1 - t**0.8)  # G: fades
                colors[flags, 2] = 0.15 + 0.70 * t**0.6  # B: rises → magenta
                colors[flags, 3] = 0.65 + 0.35 * t  # A: brighter when diving
            elif flags.any():
                colors[flags] = [0.9, 0.2, 0.5, 0.85]  # flat magenta fallback
            return colors

        if mode == "roi" and "IN_ROI" in df.columns:
            flags = df["IN_ROI"].values[idx].astype(bool)
            hi = np.array([1, 0.9, 0, 1])
            lo = np.array([0.5, 0.5, 0.5, 0.12])
            return np.where(flags[:, None], hi, lo)

        # ── Diverging mode ──
        if mode == "radial_vel" and "RADIAL_VELOCITY_SMOOTH" in df.columns:
            v = df["RADIAL_VELOCITY_SMOOTH"].values[idx].copy()
            nan_mask = np.isnan(v)
            v[nan_mask] = 0.0
            valid = v[~nan_mask]
            vmax = (
                max(np.abs(np.nanpercentile(valid, [2, 98])).max(), 1e-6)
                if len(valid) > 0
                else 1.0
            )
            v_clip = np.clip(v / vmax, -1, 1)
            colors = np.zeros((len(v), 4))
            colors[:, 3] = 1.0
            neg = v_clip < 0
            pos = v_clip >= 0
            t_neg = 1 + v_clip[neg]
            colors[neg, 0] = t_neg
            colors[neg, 1] = t_neg
            colors[neg, 2] = 1.0
            t_pos = v_clip[pos]
            colors[pos, 0] = 1.0
            colors[pos, 1] = 1.0 - t_pos
            colors[pos, 2] = 1.0 - t_pos
            colors[nan_mask] = [0.5, 0.5, 0.5, 0.3]
            return colors

        # ── Continuous modes ──
        if mode == "depth" and "SPHERICAL_DEPTH" in df.columns:
            v = df["SPHERICAL_DEPTH"].values[idx]
        elif mode == "depth_travel" and "SPHERICAL_DEPTH" in df.columns and "TRACK_ID" in df.columns:
            # Per-track total depth range: max - min SPHERICAL_DEPTH
            depth_range = df.groupby("TRACK_ID")["SPHERICAL_DEPTH"].transform(
                lambda s: s.max() - s.min()
            )
            v = depth_range.values[idx].astype(float)
            # Spots without a track get NaN; give them grey
            nan_mask = np.isnan(v)
            v[nan_mask] = 0.0
            # Percentile clipping for better contrast
            valid = v[~nan_mask]
            if len(valid) > 0:
                vmax = max(np.percentile(valid, 98), 1e-6)
                v_norm = np.clip(v / vmax, 0, 1)
            else:
                v_norm = np.zeros_like(v)
            colors = _viridis_colors(v_norm)
            colors[nan_mask] = [0.4, 0.4, 0.4, 0.15]
            return colors
        elif mode == "theta" and "THETA_DEG" in df.columns:
            v = df["THETA_DEG"].values[idx]
        elif mode == "phi" and "PHI_DEG" in df.columns:
            v = df["PHI_DEG"].values[idx]
        elif mode == "speed":
            v = self._compute_spot_speed()[idx]
        else:
            # Default: colour by frame
            v = df["FRAME"].values[idx].astype(float)

        v = v.astype(float)
        v_norm = (v - np.nanmin(v)) / max(np.nanmax(v) - np.nanmin(v), 1e-10)
        return _viridis_colors(v_norm)

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
        #
        # Background layer: all non-ingressing tracks rendered in faint
        # uniform grey so they provide spatial context without distracting.
        #
        # Foreground layer: ingressing tracks colored by ONSET FRAME
        # using a warm 'inferno' colormap.  Early ingressors appear in
        # deep red/orange while later ones glow yellow/white, creating a
        # beautiful temporal accumulation effect as you scrub through time.
        # Ghost-point extension keeps every track visible through the
        # movie's last frame.
        #
        if mode == "ingressing" and "INGRESSING" in df_t.columns:
            ing_mask = df_t["INGRESSING"].values.astype(bool)

            # ── Background: non-ingressing tracks ──
            if (~ing_mask).any():
                df_other = df_t[~ing_mask]
                if len(df_other) >= 3:
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
                            tail_width=max(0.5, self._track_width * 0.2),
                            tail_length=15,
                            color_by="track_id",
                            colormap="gray",
                            opacity=0.10,
                        )

            # ── Foreground: ingressing tracks (onset-coloured) ──
            if ing_mask.any():
                df_ing = df_t[ing_mask]
                tc = df_ing.groupby("TRACK_ID")["FRAME"].transform("count")
                df_ing = df_ing[tc >= 2].sort_values(["TRACK_ID", "FRAME"])
                if len(df_ing) > 0:
                    # Ghost-point extension: add synthetic point at the
                    # dataset's maximum frame so ingressing tracks persist
                    # visually through the end of the recording.
                    max_frame = int(self.spots["FRAME"].max())
                    last_pts = df_ing.groupby("TRACK_ID").tail(1).copy()
                    needs_ext = last_pts[last_pts["FRAME"] < max_frame]
                    if len(needs_ext) > 0:
                        ext = needs_ext.copy()
                        ext["FRAME"] = max_frame
                        df_ing = pd.concat([df_ing, ext], ignore_index=True)
                        df_ing = df_ing.sort_values(["TRACK_ID", "FRAME"])

                    t_ing = np.column_stack(
                        [
                            df_ing["TRACK_ID"].values,
                            df_ing["FRAME"].values,
                            df_ing["POSITION_Z"].values,
                            df_ing["POSITION_Y"].values,
                            df_ing["POSITION_X"].values,
                        ]
                    )
                    n_frames = (
                        int(self.spots["FRAME"].nunique())
                        if self.spots is not None
                        else 9999
                    )

                    # Color by radial velocity → emphasises the
                    # ingression movement itself.  Each point along the
                    # track is colored by how fast the cell dives inward
                    # at that moment: gentle inward = warm orange,
                    # strong dive = bright magenta/pink.
                    kw_ing = dict(
                        name="Ingressing Tracks",
                        tail_width=self._track_width,
                        tail_length=n_frames,
                    )
                    if "RADIAL_VELOCITY_SMOOTH" in df_ing.columns:
                        vr = df_ing["RADIAL_VELOCITY_SMOOTH"].values.copy()
                        vr = np.nan_to_num(vr, nan=0.0)
                        # Inward speed: flip sign so positive = diving in
                        inward = np.clip(-vr, 0, None)
                        vmax = (
                            np.percentile(inward[inward > 0], 95)
                            if (inward > 0).any()
                            else 1.0
                        )
                        # Normalize 0–1 for colormap
                        normed = np.clip(inward / max(vmax, 1e-6), 0, 1)
                        kw_ing["properties"] = {"inward_speed": normed}
                        kw_ing["color_by"] = "inward_speed"
                        kw_ing["colormap"] = "magma"
                    else:
                        kw_ing["color_by"] = "track_id"
                        kw_ing["colormap"] = "magenta"

                    self.tracks_layer = self.viewer.add_tracks(t_ing, **kw_ing)
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
        elif mode == "depth_travel" and "SPHERICAL_DEPTH" in df_t.columns:
            # Per-track total depth range: colour each track by how far
            # it travels radially through the blastoderm
            depth_range = df_t.groupby("TRACK_ID")["SPHERICAL_DEPTH"].transform(
                lambda s: s.max() - s.min()
            ).values.astype(float)
            nan_m = np.isnan(depth_range)
            if not nan_m.all():
                depth_range[nan_m] = 0.0
                # Percentile clipping for better contrast across all tracks
                valid = depth_range[~nan_m]
                if len(valid) > 0:
                    vmax = max(np.percentile(valid, 98), 1e-6)
                    depth_range = np.clip(depth_range / vmax, 0, 1)
                props["color_value"] = depth_range
                color_by = "color_value"
                cmap = "plasma"
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
            tail_width=self._track_width,
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

            # Update max tracks slider to actual track count
            if "TRACK_ID" in self.spots.columns:
                n_tracks = max(500, int(self.spots["TRACK_ID"].nunique()))
                self.max_tracks_slider.max = n_tracks
                self.max_tracks_slider.value = n_tracks
                self._max_tracks = n_tracks

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
            fmt = _detect_csv_format(tracks_path)
            if fmt == "ultrack":
                self.tracks_raw = load_ultrack_csv(tracks_path)
            else:
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

        btn_load_segments = PushButton(text="Load Segments (.zarr)")
        btn_load_segments.changed.connect(self._load_segments_from_ui)
        self.lbl_segments_status = Label(
            value=f"Loaded: {self._segments_data.shape}"
            if self._segments_data is not None
            else "No segments loaded"
        )
        self.segments_visible_cb = CheckBox(value=True, text="Show 3D nuclei")
        self.segments_visible_cb.changed.connect(self._on_segments_visible_changed)
        self.segments_opacity_slider = FloatSlider(
            value=50, min=0, max=100, step=5, label="Nuclei 3D opacity %"
        )
        self.segments_opacity_slider.changed.connect(self._on_segments_opacity_changed)

        # ── Display controls ──
        self.display_pct_slider = FloatSlider(
            value=100, min=1, max=100, step=1, label="Display %"
        )
        self.display_pct_slider.changed.connect(self._on_display_pct_changed)

        n_tracks_total = 500000
        if self.spots is not None and "TRACK_ID" in self.spots.columns:
            n_tracks_total = max(500, int(self.spots["TRACK_ID"].nunique()))
        self.max_tracks_slider = FloatSlider(
            value=n_tracks_total,
            min=100,
            max=n_tracks_total,
            step=500,
            label="Max tracks",
        )
        self._max_tracks = n_tracks_total
        self.max_tracks_slider.changed.connect(self._on_max_tracks_changed)

        self.track_width_slider = FloatSlider(
            value=2.0, min=0.5, max=10.0, step=0.5, label="Track width"
        )
        self.track_width_slider.changed.connect(self._on_track_width_changed)

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

        btn_margin = PushButton(text="Pick Margin")
        btn_margin.changed.connect(lambda: self._set_pick_mode("margin"))

        btn_ingression = PushButton(text="Pick Ingression Center")
        btn_ingression.changed.connect(lambda: self._set_pick_mode("ingression"))

        self.lbl_ap = Label(value="AP: not set")
        self.lbl_dorsal = Label(value="Dorsal: not set")
        self.lbl_margin = Label(value="Margin: not set")
        self.lbl_ingression = Label(value="Ingression: not set")

        # ── Orient ──
        btn_orient = PushButton(text="Orient Embryo")
        btn_orient.changed.connect(self._orient)

        # ── Sphere fit ──
        btn_sphere = PushButton(text="Fit Sphere")
        btn_sphere.changed.connect(self._fit_sphere)
        btn_pick_sphere_center = PushButton(text="📍 Pick Sphere Center")
        btn_pick_sphere_center.changed.connect(
            lambda: self._set_pick_mode("sphere_center")
        )
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

        self.ingression_roi_only = CheckBox(value=False, text="Flag only within ROI")

        # ── Colour mode ──
        btn_color_frame = PushButton(text="Colour by Frame")
        btn_color_frame.changed.connect(lambda: self._recolor("frame"))
        btn_color_depth = PushButton(text="Colour by Depth")
        btn_color_depth.changed.connect(lambda: self._recolor("depth"))
        btn_color_depth_travel = PushButton(text="Colour by Depth Travel")
        btn_color_depth_travel.changed.connect(lambda: self._recolor("depth_travel"))
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

        # ── Cross-section ──
        btn_cross_section = PushButton(text="✂ Cross-Section View")
        btn_cross_section.changed.connect(self._open_cross_section)

        container = Container(
            widgets=[
                Label(value="── Data ──"),
                btn_load_spots,
                self.lbl_spots_status,
                btn_load_tracks,
                self.lbl_tracks_status,
                btn_load_segments,
                self.lbl_segments_status,
                self.segments_visible_cb,
                self.segments_opacity_slider,
                Label(value="━━ ⏱ TIME WINDOW ━━"),
                self.track_frame_start,
                self.track_frame_end,
                Label(value="── Display ──"),
                self.display_pct_slider,
                self.max_tracks_slider,
                self.track_width_slider,
                Label(value="── Landmarks ──"),
                btn_ap,
                self.lbl_ap,
                btn_dorsal,
                self.lbl_dorsal,
                btn_margin,
                self.lbl_margin,
                btn_ingression,
                self.lbl_ingression,
                Label(value="── Orientation ──"),
                btn_orient,
                Label(value="── Sphere Fit ──"),
                btn_sphere,
                btn_pick_sphere_center,
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
                self.ingression_roi_only,
                btn_flag_ingression,
                self.lbl_ingression,
                self.ingression_tracks_only,
                Label(value="── Colouring ──"),
                btn_color_frame,
                btn_color_depth,
                btn_color_depth_travel,
                btn_color_theta,
                btn_color_phi,
                btn_color_speed,
                btn_color_radvel,
                btn_color_ingress,
                btn_color_roi,
                Label(value="── Export & Video ──"),
                btn_export,
                btn_record_video,
                btn_cross_section,
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
                "sphere_center": "SPHERE CENTER",
                "margin": "MARGIN (click a cell at the margin boundary)",
                "ingression": "INGRESSION CENTER (click where ingression is strongest)",
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
        elif self._pick_mode == "sphere_center":
            self._pick_mode = "none"
            self._apply_manual_sphere_center(pos_xyz)
        elif self._pick_mode == "margin":
            self.margin_landmark = pos_xyz
            # Compute θ if sphere is fitted
            theta_str = ""
            if self.sphere_fit_result is not None:
                c = self.sphere_fit_result["center"]
                d = pos_xyz - c
                r = np.linalg.norm(d)
                theta = np.degrees(np.arccos(np.clip(-d[1] / max(r, 1e-10), -1, 1)))
                theta_str = f" | θ={theta:.1f}°"
            self.lbl_margin.value = (
                f"Margin: Y={pos_xyz[1]:.1f}{theta_str}"
            )
            self._pick_mode = "none"
            self.viewer.status = (
                f"Margin set at Y={pos_xyz[1]:.1f}{theta_str} "
                f"(snapped {snap_d:.1f}µm to nearest spot)"
            )
        elif self._pick_mode == "ingression":
            self.ingression_landmark = pos_xyz
            theta_str, phi_str = "", ""
            if self.sphere_fit_result is not None:
                c = self.sphere_fit_result["center"]
                d = pos_xyz - c
                r = np.linalg.norm(d)
                theta = np.degrees(np.arccos(np.clip(-d[1] / max(r, 1e-10), -1, 1)))
                phi = np.degrees(np.arctan2(d[2], d[0]))
                theta_str = f" | θ={theta:.1f}°"
                phi_str = f", φ={phi:.1f}°"
            self.lbl_ingression.value = (
                f"Ingression: ({pos_xyz[0]:.1f}, {pos_xyz[1]:.1f}, {pos_xyz[2]:.1f})"
                f"{theta_str}{phi_str}"
            )
            self._pick_mode = "none"
            self.viewer.status = (
                f"Ingression center set{theta_str}{phi_str} "
                f"(snapped {snap_d:.1f}µm to nearest spot)"
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
        if self.margin_landmark is not None:
            pts.append([self.margin_landmark[2], self.margin_landmark[1], self.margin_landmark[0]])
            colors.append("yellow")
        if self.ingression_landmark is not None:
            pts.append([self.ingression_landmark[2], self.ingression_landmark[1], self.ingression_landmark[0]])
            colors.append("magenta")

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

    def _apply_manual_sphere_center(self, center_xyz):
        """Set sphere center manually; estimate radius from data."""
        x = self.spots["POSITION_X"].values
        y = self.spots["POSITION_Y"].values
        z = self.spots["POSITION_Z"].values

        pts = np.column_stack([x, y, z])
        dists = np.sqrt(np.sum((pts - center_xyz) ** 2, axis=1))
        radius = float(np.median(dists))
        rmse = float(np.sqrt(np.mean((dists - radius) ** 2)))

        result = {
            "center": center_xyz,
            "radius": radius,
            "rmse": rmse,
            "coverage": 1.0,
            "surface_used": len(pts),
            "cap_height": radius,
            "cap_base_radius": radius,
        }
        self._finish_sphere_fit(result)
        self.viewer.status = (
            f"Manual sphere center set — R={radius:.0f}µm, RMSE={rmse:.1f}µm"
        )

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

        # Validate mode is available
        df = self.spots
        needed = {
            "depth": "SPHERICAL_DEPTH",
            "depth_travel": "SPHERICAL_DEPTH",
            "theta": "THETA_DEG",
            "phi": "PHI_DEG",
            "radial_vel": "RADIAL_VELOCITY_SMOOTH",
            "ingressing": "INGRESSING",
            "roi": "IN_ROI",
        }
        if mode in needed and needed[mode] not in df.columns:
            self.viewer.status = (
                f"'{mode}' not available — fit sphere / apply ROI first."
            )
            return

        self._current_color_mode = mode
        colors = self._get_display_colors()
        self.spots_layer.face_color = colors
        self._rebuild_tracks()
        self.viewer.status = f"Coloured by {mode}"

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

        # If "Flag only within ROI" is checked, apply two spatial filters:
        #
        #  (1) ROI PRESENCE — the track must have ≥1 spot inside the ROI
        #      box.  This removes false positives near the animal pole
        #      (entirely above the box) while preserving tracks that start
        #      outside the ROI and later migrate down into the marginal
        #      zone during ingression.
        #
        #  (2) VEGETAL BREACH — any track where ANY spot crosses past the
        #      vegetal boundary (y_max, after orientation +Y = vegetal) is
        #      epiboly / spreading, not ingression, so it is excluded.
        #
        if self.ingression_roi_only.value and self._roi_bounds is not None:
            b = self._roi_bounds
            df = self.spots

            # (1) Track must touch the ROI at least once
            in_box = (
                (df["POSITION_X"] >= b["x_min"])
                & (df["POSITION_X"] <= b["x_max"])
                & (df["POSITION_Y"] >= b["y_min"])
                & (df["POSITION_Y"] <= b["y_max"])
            )
            track_ever_in_box = in_box.groupby(df["TRACK_ID"]).transform("any")
            self.spots.loc[~track_ever_in_box, "INGRESSING"] = False

            # (2) Exclude tracks breaching the vegetal boundary
            vegetal_breach = (
                df.groupby("TRACK_ID")["POSITION_Y"].transform("max") > b["y_max"]
            )
            self.spots.loc[vegetal_breach, "INGRESSING"] = False

        n_ing = int(self.spots["INGRESSING"].sum())
        n_total = len(self.spots)
        n_tracks = 0
        if "TRACK_ID" in self.spots.columns:
            n_tracks = int(
                self.spots.loc[self.spots["INGRESSING"], "TRACK_ID"].nunique()
            )
        pct = n_ing / max(n_total, 1) * 100
        roi_tag = " [ROI-only]" if self.ingression_roi_only.value else ""
        self.lbl_ingression.value = (
            f"Ingressing: {n_ing:,} spots ({n_tracks} tracks, {pct:.1f}%){roi_tag}\n"
            f"V_r≤{thresh:.3f}, inward≥{inward_frac:.0%}, ≥{min_f}f, net-in+sustained"
        )
        self.viewer.status = (
            f"Flagged {n_ing:,} ingressing spots ({n_tracks} tracks){roi_tag}. "
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
        data, idx = self._build_display_points()
        self._display_idx = idx
        colors = self._get_display_colors()
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

    def _on_track_width_changed(self):
        """Handle track width slider change (debounced)."""
        self._track_width = float(self.track_width_slider.value)
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

    # ── Cross-section view ────────────────────────────────────────────

    def _open_cross_section(self):
        """Open a cross-section 2D scatter viewer in a separate window."""
        if self.spots is None:
            self.viewer.status = "Load spots first!"
            return

        self.viewer.status = "Opening cross-section viewer..."
        CrossSectionViewer(self)

    # ── Export & Video ────────────────────────────────────────────────

    def _export(self):
        """Export comprehensive analysis output for downstream R pipeline.

        Runs the heavy CSV writes in a background thread to avoid freezing
        the UI on large datasets (millions of rows).
        """
        import json

        if self.spots is None:
            self.viewer.status = "Load spots first!"
            return

        from qtpy.QtCore import QTimer

        out_dir = Path("analysis_output")
        out_dir.mkdir(exist_ok=True)

        self.lbl_export.value = "Exporting (please wait)..."
        self.viewer.status = "Exporting analysis files in background..."

        # Snapshot everything we need — avoids touching GUI state from thread
        spots_df = self.spots.copy()
        tracks_raw = self.tracks_raw.copy() if self.tracks_raw is not None else None
        sphere = (
            dict(self.sphere_fit_result) if self.sphere_fit_result is not None else None
        )
        roi = dict(self._roi_bounds) if self._roi_bounds is not None else None
        roi_flag = getattr(self, "ingression_roi_only", None)
        roi_restricted = roi_flag is not None and roi_flag.value
        thresh_val = self.ingression_threshold.value
        min_f_val = int(self.ingression_min_frames.value)
        inward_frac_val = self.ingression_inward_frac.value
        oriented = self.oriented
        track_frame_min = self._track_frame_min
        track_frame_max = self._track_frame_max
        display_pct = self._display_pct
        max_tracks = self._max_tracks
        color_mode = self._current_color_mode

        # Snapshot gastrulation landmarks (margin + ingression)
        margin_lm = self.margin_landmark.copy() if self.margin_landmark is not None else None
        ingression_lm = self.ingression_landmark.copy() if self.ingression_landmark is not None else None
        sphere_for_lm = sphere  # already snapshot

        self._pending_export_err = None
        self._pending_export_files = None

        def _bg():
            try:
                files_written = []

                # 1. Enriched spots CSV (all computed columns)
                spots_df.to_csv(out_dir / "oriented_spots.csv", index=False)
                files_written.append("oriented_spots.csv")

                # 2. Enriched tracks CSV — spots filtered to tracked nuclei
                #    only, sorted by TRACK_ID + FRAME for easy downstream use
                if "TRACK_ID" in spots_df.columns:
                    tracked = spots_df.dropna(subset=["TRACK_ID"]).sort_values(
                        ["TRACK_ID", "FRAME"]
                    )
                    if len(tracked) > 0:
                        tracked.to_csv(out_dir / "oriented_tracks.csv", index=False)
                        files_written.append("oriented_tracks.csv")

                # 3. Sphere parameters
                if sphere is not None:
                    c = sphere["center"]
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
                                sphere["radius"],
                                sphere["rmse"],
                                sphere["coverage"],
                                sphere["surface_used"],
                                sphere["cap_height"],
                                sphere["cap_base_radius"],
                            ],
                        }
                    )
                    params.to_csv(out_dir / "sphere_params.csv", index=False)
                    files_written.append("sphere_params.csv")

                # 4. ROI bounds
                if roi is not None:
                    roi_df = pd.DataFrame(
                        {"parameter": list(roi.keys()), "value": list(roi.values())}
                    )
                    roi_df.to_csv(out_dir / "roi_bounds.csv", index=False)
                    files_written.append("roi_bounds.csv")

                # 5. Ingression parameters
                if "INGRESSING" in spots_df.columns:
                    n_ing = int(spots_df["INGRESSING"].sum())
                    n_tracks_ing = 0
                    if "TRACK_ID" in spots_df.columns:
                        n_tracks_ing = int(
                            spots_df.loc[
                                spots_df["INGRESSING"].astype(bool), "TRACK_ID"
                            ].nunique()
                        )
                    ing_df = pd.DataFrame(
                        {
                            "parameter": [
                                "threshold_um_per_frame",
                                "min_track_frames",
                                "min_inward_fraction",
                                "roi_restricted",
                                "detection_method",
                                "criteria",
                                "n_ingressing_spots",
                                "n_ingressing_tracks",
                            ],
                            "value": [
                                thresh_val,
                                min_f_val,
                                inward_frac_val,
                                roi_restricted,
                                "embryology-informed multi-criteria",
                                "median_Vr+inward_frac+net_disp+sustained",
                                n_ing,
                                n_tracks_ing,
                            ],
                        }
                    )
                    ing_df.to_csv(out_dir / "ingression_params.csv", index=False)
                    files_written.append("ingression_params.csv")

                # 5b. Gastrulation landmarks (margin + ingression center)
                lm_rows = []
                if margin_lm is not None:
                    lm_rows.append({"landmark": "margin", "x": margin_lm[0],
                                    "y": margin_lm[1], "z": margin_lm[2]})
                    if sphere_for_lm is not None:
                        c = sphere_for_lm["center"]
                        d = margin_lm - c
                        r = np.linalg.norm(d)
                        theta = float(np.degrees(np.arccos(np.clip(-d[1] / max(r, 1e-10), -1, 1))))
                        phi = float(np.degrees(np.arctan2(d[2], d[0])))
                        lm_rows[-1].update({"theta": theta, "phi": phi})
                if ingression_lm is not None:
                    lm_rows.append({"landmark": "ingression_center", "x": ingression_lm[0],
                                    "y": ingression_lm[1], "z": ingression_lm[2]})
                    if sphere_for_lm is not None:
                        c = sphere_for_lm["center"]
                        d = ingression_lm - c
                        r = np.linalg.norm(d)
                        theta = float(np.degrees(np.arccos(np.clip(-d[1] / max(r, 1e-10), -1, 1))))
                        phi = float(np.degrees(np.arctan2(d[2], d[0])))
                        lm_rows[-1].update({"theta": theta, "phi": phi})
                if lm_rows:
                    lm_df = pd.DataFrame(lm_rows)
                    lm_df.to_csv(out_dir / "gastrulation_landmarks.csv", index=False)
                    files_written.append("gastrulation_landmarks.csv")

                # 6. Per-track summary CSV
                if "TRACK_ID" in spots_df.columns:
                    track_cols = [
                        "TRACK_MEDIAN_RADIAL_VEL",
                        "TRACK_MEAN_RADIAL_VEL",
                        "TRACK_INWARD_FRACTION",
                        "TRACK_MAX_INWARD_VEL",
                        "TRACK_NET_RADIAL_DISP",
                        "TRACK_SUSTAINED_INWARD",
                        "INGRESSION_ONSET_FRAME",
                        "INGRESSING",
                        "INGRESSION_SCORE",
                    ]
                    avail_cols = [c for c in track_cols if c in spots_df.columns]
                    if avail_cols:
                        grp = spots_df.dropna(subset=["TRACK_ID"]).groupby("TRACK_ID")
                        track_summary = grp.agg(
                            n_spots=("FRAME", "count"),
                            frame_start=("FRAME", "min"),
                            frame_end=("FRAME", "max"),
                            mean_x=("POSITION_X", "mean"),
                            mean_y=("POSITION_Y", "mean"),
                            mean_z=("POSITION_Z", "mean"),
                        )
                        first_per_track = grp[avail_cols].first()
                        track_summary = track_summary.join(first_per_track)
                        for sc in ["SPHERICAL_DEPTH", "THETA_DEG", "PHI_DEG"]:
                            if sc in spots_df.columns:
                                track_summary[f"mean_{sc}"] = grp[sc].mean()
                        if "IN_ROI" in spots_df.columns:
                            track_summary["IN_ROI"] = grp["IN_ROI"].any()
                        track_summary.to_csv(out_dir / "track_summary.csv")
                        files_written.append("track_summary.csv")

                # 7. Analysis metadata (JSON)
                meta = {
                    "oriented": oriented,
                    "sphere_fitted": sphere is not None,
                    "roi_set": roi is not None,
                    "roi_restricted_ingression": roi_restricted,
                    "time_window": {
                        "frame_min": track_frame_min,
                        "frame_max": track_frame_max,
                    },
                    "display": {
                        "display_pct": display_pct,
                        "max_tracks": max_tracks,
                        "current_color_mode": color_mode,
                    },
                    "n_spots": len(spots_df),
                    "n_frames": int(spots_df["FRAME"].nunique()),
                    "columns_exported": list(spots_df.columns),
                    "ingression_criteria": {
                        "method": "multi-criteria embryology-informed",
                        "criteria": [
                            "median V_r <= threshold",
                            "inward fraction >= min_inward_fraction",
                            "track length >= min_track_frames",
                            "net radial displacement >= 0.5 µm inward",
                            "sustained consecutive inward >= max(3, 15% of track)",
                        ],
                    },
                }
                if sphere is not None:
                    meta["sphere_radius"] = float(sphere["radius"])
                if roi is not None:
                    meta["roi_bounds"] = {k: float(v) for k, v in roi.items()}
                with open(out_dir / "analysis_metadata.json", "w") as f:
                    json.dump(meta, f, indent=2, default=str)
                files_written.append("analysis_metadata.json")

                self._pending_export_files = files_written
            except Exception as e:
                self._pending_export_err = str(e)

        t = threading.Thread(target=_bg, daemon=True)
        t.start()

        def _poll():
            if t.is_alive():
                QTimer.singleShot(500, _poll)
                return
            if self._pending_export_err:
                self.lbl_export.value = f"Export error: {self._pending_export_err}"
                self.viewer.status = f"Export failed: {self._pending_export_err}"
                return
            files = self._pending_export_files or []
            self.lbl_export.value = f"Exported {len(files)} files to {out_dir}/"
            self.viewer.status = f"Exported: {', '.join(files)}"

        QTimer.singleShot(500, _poll)

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
                # Convert RGBA → RGB and ensure even dimensions for libx264
                if img.ndim == 3 and img.shape[2] == 4:
                    img = img[:, :, :3]
                h, w = img.shape[:2]
                if h % 2 or w % 2:
                    img = img[: h - h % 2, : w - w % 2]
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
        description="napari-based 4D embryo viewer (TrackMate & ultrack formats)"
    )
    parser.add_argument(
        "spots",
        nargs="?",
        default=None,
        help="Path to spots/tracks CSV — auto-detects TrackMate or ultrack format",
    )
    parser.add_argument(
        "--tracks", "-t", default=None, help="Path to tracks CSV (optional)"
    )
    parser.add_argument(
        "--segments", "-s", default=None,
        help="Path to segments.zarr (ultrack label volume, 4D)",
    )
    parser.add_argument(
        "--voxel-size", type=float, nargs=3, default=[1.0, 1.0, 1.0],
        metavar=("Z", "Y", "X"),
        help="Voxel size in µm (Z Y X) for segments layer. Default: 1.0 1.0 1.0",
    )
    parser.add_argument(
        "--max-spots",
        type=int,
        default=None,
        help="Max spots to load (subsample for speed)",
    )
    parser.add_argument(
        "--downsample",
        type=int,
        default=2,
        metavar="N",
        help="Spatial downsample factor for segments display layer "
             "(e.g. 2 = half-res → 8× fewer voxels, 4 = quarter-res → "
             "64× fewer). Full-res data is kept in RAM for analysis. "
             "Default: 2 (good balance of quality and scrolling speed). "
             "Use 1 for full-res display if GPU can handle it.",
    )
    parser.add_argument(
        "--preload",
        action="store_true",
        help="Load entire segments.zarr into RAM with uint16 remapping "
             "and GPU-optimised display volume for instant smooth "
             "scrubbing. Recommended for machines with enough RAM.",
    )
    args = parser.parse_args()

    spots_df = None
    tracks_df = None

    # Determine spots source: explicit spots arg, or fall back to --tracks
    spots_path = args.spots
    if spots_path is None and args.tracks is not None:
        spots_path = args.tracks  # tracks CSV has all position columns

    if spots_path:
        fmt = _detect_csv_format(spots_path)
        print(f"Loading spots from {spots_path} (detected: {fmt} format)...")
        spots_df = load_spots(spots_path)
        print(f"  {len(spots_df):,} spots, {spots_df['FRAME'].nunique()} frames")

        if args.max_spots and len(spots_df) > args.max_spots:
            print(f"  Subsampling to {args.max_spots:,} spots...")
            spots_df = spots_df.sample(args.max_spots, random_state=42)

    if args.tracks:
        print(f"Loading tracks from {args.tracks}...")
        fmt_t = _detect_csv_format(args.tracks)
        if fmt_t == "ultrack":
            tracks_df = load_ultrack_csv(args.tracks)
        else:
            tracks_df = load_trackmate_csv(args.tracks)
        print(f"  {len(tracks_df):,} track records")

    if spots_df is None and args.segments is None:
        print("Launching viewer — use buttons to open files.")

    print("Launching napari viewer...")
    viewer = EmbryoViewer(
        spots_df, tracks_df,
        segments_zarr_path=args.segments,
        voxel_size=tuple(args.voxel_size),
        downsample=args.downsample,
        preload=args.preload,
    )

    napari, _ = _import_napari()
    napari.run()


if __name__ == "__main__":
    main()
