#!/usr/bin/env python3
"""
Compute nuclear statistics from segmented label images (per-timepoint TIFs).

For each timepoint and z-slice, measures:
  - Nuclear density   (nuclei / area in µm²)
  - Nuclear size      (area in µm² per nucleus)
  - Internuclear dist (nearest-neighbour distance in µm between centroids)

Also produces per-timepoint 3D summary (volume-based density, 3D NN distance).

Outputs:
  nuclear_stats_per_slice.tsv   — per (timepoint, z-slice) stats
  nuclear_stats_per_time.tsv    — per timepoint (3D) stats
  nuclear_stats_global.tsv      — global summary across all timepoints

Usage (on the cluster):
  python nuclear_stats.py /path/to/segmented -o /path/to/output

The segmented folder should contain per-timepoint TIFs named like
t0001_segmented.tif, t0002_segmented.tif, …  (3D label stacks, ZxYxX).

Voxel size is auto-read from TIFF metadata when available; otherwise
use --voxel-size Z Y X (in µm).
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import label as ndi_label
from scipy.spatial import KDTree
from skimage.measure import regionprops_table
from tifffile import TiffFile


# ──────────────────────────────────────────────────────────────────────────────
#  Helpers
# ──────────────────────────────────────────────────────────────────────────────

def _natural_sort_key(p: Path) -> tuple:
    """Sort key that handles embedded integers (t0001 < t0010 < t0100)."""
    parts = re.split(r"(\d+)", p.stem)
    return tuple(int(x) if x.isdigit() else x.lower() for x in parts)


def read_voxel_size(tif_path: Path) -> tuple[float, float, float]:
    """Try to extract (dz, dy, dx) in µm from TIFF metadata.

    Falls back to (1.0, 1.0, 1.0) if nothing useful is found.
    """
    with TiffFile(str(tif_path)) as tif:
        # ImageJ hyperstack metadata
        if tif.imagej_metadata and "spacing" in tif.imagej_metadata:
            dz = float(tif.imagej_metadata["spacing"])
        else:
            dz = 1.0

        # XY resolution from first page
        page = tif.pages[0]
        tags = page.tags
        if "XResolution" in tags and "YResolution" in tags:
            xr = tags["XResolution"].value
            yr = tags["YResolution"].value
            # TIFF stores resolution as (numerator, denominator) — pixels/unit
            dx = xr[1] / xr[0] if xr[0] != 0 else 1.0
            dy = yr[1] / yr[0] if yr[0] != 0 else 1.0
            # ResolutionUnit: 1=no unit, 2=inch, 3=cm
            unit_tag = tags.get("ResolutionUnit")
            unit = unit_tag.value if unit_tag else 1
            if unit == 2:  # inch → µm
                dx *= 25_400
                dy *= 25_400
            elif unit == 3:  # cm → µm
                dx *= 10_000
                dy *= 10_000
        else:
            dx, dy = 1.0, 1.0

    return (dz, dy, dx)


def _extract_timepoint(path: Path) -> int:
    """Extract integer timepoint from filename like t0001_segmented.tif."""
    m = re.search(r"t(\d+)", path.stem, re.IGNORECASE)
    if m:
        return int(m.group(1))
    # fallback: first integer found
    m = re.search(r"(\d+)", path.stem)
    return int(m.group(1)) if m else 0


# ──────────────────────────────────────────────────────────────────────────────
#  Per-slice statistics
# ──────────────────────────────────────────────────────────────────────────────

def compute_slice_stats(
    label_slice: np.ndarray,
    pixel_area_um2: float,
) -> dict:
    """Compute stats for a single 2D label slice.

    Parameters
    ----------
    label_slice : 2D ndarray of integer labels (0 = background).
    pixel_area_um2 : area of one pixel in µm².

    Returns
    -------
    dict with keys: n_nuclei, slice_area_um2, density_per_um2,
                    size_mean_um2, size_median_um2, size_std_um2,
                    nn_dist_mean_um, nn_dist_median_um, nn_dist_std_um
    """
    unique_labels = np.unique(label_slice)
    unique_labels = unique_labels[unique_labels != 0]
    n_nuclei = len(unique_labels)

    slice_area_um2 = label_slice.shape[0] * label_slice.shape[1] * pixel_area_um2

    result = {
        "n_nuclei": n_nuclei,
        "slice_area_um2": slice_area_um2,
        "density_per_um2": n_nuclei / slice_area_um2 if slice_area_um2 > 0 else np.nan,
    }

    if n_nuclei == 0:
        result.update({
            "size_mean_um2": np.nan,
            "size_median_um2": np.nan,
            "size_std_um2": np.nan,
            "eqcircle_diam_mean_um": np.nan,
            "eqcircle_diam_median_um": np.nan,
            "eqcircle_diam_std_um": np.nan,
            "nn_dist_mean_um": np.nan,
            "nn_dist_median_um": np.nan,
            "nn_dist_std_um": np.nan,
        })
        return result

    # Region properties (area + centroid)
    props = regionprops_table(
        label_slice,
        properties=("label", "area", "centroid"),
    )
    areas_um2 = np.array(props["area"]) * pixel_area_um2
    result["size_mean_um2"] = float(np.mean(areas_um2))
    result["size_median_um2"] = float(np.median(areas_um2))
    result["size_std_um2"] = float(np.std(areas_um2))

    # Equivalent-circle diameter: d = 2 * sqrt(A / π)
    eqcircle_diams = 2.0 * np.sqrt(areas_um2 / np.pi)
    result["eqcircle_diam_mean_um"] = float(np.mean(eqcircle_diams))
    result["eqcircle_diam_median_um"] = float(np.median(eqcircle_diams))
    result["eqcircle_diam_std_um"] = float(np.std(eqcircle_diams))

    # Nearest-neighbour distances (in µm, using pixel_size for Y, X)
    centroids = np.column_stack([props["centroid-0"], props["centroid-1"]])
    if n_nuclei >= 2:
        # Scale centroids to µm (assume square pixels in 2D — use sqrt(pixel_area))
        px_size = np.sqrt(pixel_area_um2)
        centroids_um = centroids * px_size
        tree = KDTree(centroids_um)
        dists, _ = tree.query(centroids_um, k=2)
        nn_dists = dists[:, 1]  # nearest neighbour (skip self)
        result["nn_dist_mean_um"] = float(np.mean(nn_dists))
        result["nn_dist_median_um"] = float(np.median(nn_dists))
        result["nn_dist_std_um"] = float(np.std(nn_dists))
    else:
        result["nn_dist_mean_um"] = np.nan
        result["nn_dist_median_um"] = np.nan
        result["nn_dist_std_um"] = np.nan

    return result


# ──────────────────────────────────────────────────────────────────────────────
#  Per-timepoint 3D statistics
# ──────────────────────────────────────────────────────────────────────────────

CORE_Z_FRAC = 0.25  # fraction of peak nuclei count defining the core


def compute_volume_stats(
    labels_3d: np.ndarray,
    voxel_size: tuple[float, float, float],
    core_z_range: tuple[int, int] | None = None,
) -> dict:
    """Compute 3D nuclear stats for a whole stack.

    Parameters
    ----------
    labels_3d : ZxYxX label array.
    voxel_size : (dz, dy, dx) in µm.
    core_z_range : (z_lo, z_hi) inclusive slice indices defining the core.
                   If provided, core-filtered 3D stats are also computed.

    Returns
    -------
    dict with 3D summary values (all-nuclei + core-filtered if available).
    """
    dz, dy, dx = voxel_size
    voxel_vol = dz * dy * dx

    unique_labels = np.unique(labels_3d)
    unique_labels = unique_labels[unique_labels != 0]
    n_nuclei = len(unique_labels)

    vol_um3 = np.prod(labels_3d.shape) * voxel_vol

    result = {
        "n_nuclei_3d": n_nuclei,
        "volume_um3": vol_um3,
        "density_per_um3": n_nuclei / vol_um3 if vol_um3 > 0 else np.nan,
    }

    _ALL_NAN_KEYS = [
        "vol_mean_um3", "vol_median_um3", "vol_std_um3",
        "eqsph_diam_mean_um", "eqsph_diam_median_um", "eqsph_diam_std_um",
        "bbox_major_mean_um", "bbox_major_median_um", "bbox_major_std_um",
        "bbox_minor_mean_um", "bbox_minor_median_um", "bbox_minor_std_um",
        "nn3d_mean_um", "nn3d_median_um", "nn3d_std_um",
    ]
    _CORE_NAN_KEYS = [
        "core_n_nuclei_3d",
        "core_vol_mean_um3", "core_vol_median_um3",
        "core_eqsph_diam_mean_um", "core_eqsph_diam_median_um",
        "core_bbox_major_mean_um", "core_bbox_major_median_um",
        "core_bbox_minor_mean_um", "core_bbox_minor_median_um",
        "core_nn3d_mean_um", "core_nn3d_median_um",
    ]

    if n_nuclei == 0:
        for k in _ALL_NAN_KEYS + _CORE_NAN_KEYS:
            result[k] = np.nan
        return result

    # 3D region props — only cheap properties (no inertia tensor)
    props = regionprops_table(
        labels_3d,
        properties=("label", "area", "centroid", "bbox"),
    )
    # "area" in 3D = voxel count; convert to µm³
    volumes_um3 = np.array(props["area"]) * voxel_vol
    result["vol_mean_um3"] = float(np.mean(volumes_um3))
    result["vol_median_um3"] = float(np.median(volumes_um3))
    result["vol_std_um3"] = float(np.std(volumes_um3))

    # Equivalent sphere diameter  d = 2 * (3V / 4π)^(1/3)
    eqsph_diams = 2.0 * (3.0 * volumes_um3 / (4.0 * np.pi)) ** (1.0 / 3.0)
    result["eqsph_diam_mean_um"] = float(np.mean(eqsph_diams))
    result["eqsph_diam_median_um"] = float(np.median(eqsph_diams))
    result["eqsph_diam_std_um"] = float(np.std(eqsph_diams))

    # Bounding-box extents in µm (fast direct measurement)
    # bbox columns: bbox-0..2 = min(z,y,x), bbox-3..5 = max(z,y,x)
    ext_z = (np.array(props["bbox-3"]) - np.array(props["bbox-0"])) * dz
    ext_y = (np.array(props["bbox-4"]) - np.array(props["bbox-1"])) * dy
    ext_x = (np.array(props["bbox-5"]) - np.array(props["bbox-2"])) * dx
    extents = np.column_stack([ext_z, ext_y, ext_x])
    bbox_major = np.max(extents, axis=1)  # longest bbox dimension
    bbox_minor = np.min(extents, axis=1)  # shortest bbox dimension
    result["bbox_major_mean_um"] = float(np.mean(bbox_major))
    result["bbox_major_median_um"] = float(np.median(bbox_major))
    result["bbox_major_std_um"] = float(np.std(bbox_major))
    result["bbox_minor_mean_um"] = float(np.mean(bbox_minor))
    result["bbox_minor_median_um"] = float(np.median(bbox_minor))
    result["bbox_minor_std_um"] = float(np.std(bbox_minor))

    # 3D NN distance — scale centroids to µm
    centroids = np.column_stack([
        props["centroid-0"] * dz,
        props["centroid-1"] * dy,
        props["centroid-2"] * dx,
    ])
    if n_nuclei >= 2:
        tree = KDTree(centroids)
        dists, _ = tree.query(centroids, k=2)
        nn = dists[:, 1]
        result["nn3d_mean_um"] = float(np.mean(nn))
        result["nn3d_median_um"] = float(np.median(nn))
        result["nn3d_std_um"] = float(np.std(nn))
    else:
        result["nn3d_mean_um"] = np.nan
        result["nn3d_median_um"] = np.nan
        result["nn3d_std_um"] = np.nan

    # ── Core-filtered 3D stats ──────────────────────────────────────────────
    # Filter labels whose centroid-z falls within the core z-range and
    # recompute volume, equiv-sphere diameter, bbox axes, and NN distance.
    centroid_z_slices = np.array(props["centroid-0"])  # in slice coords
    if core_z_range is not None:
        core_mask = (
            (centroid_z_slices >= core_z_range[0])
            & (centroid_z_slices <= core_z_range[1])
        )
    else:
        # No core range provided — treat all as core
        core_mask = np.ones(len(centroid_z_slices), dtype=bool)

    n_core = int(core_mask.sum())
    result["core_n_nuclei_3d"] = n_core

    if n_core == 0:
        for k in _CORE_NAN_KEYS:
            if k not in result:
                result[k] = np.nan
    else:
        core_vols = volumes_um3[core_mask]
        result["core_vol_mean_um3"] = float(np.mean(core_vols))
        result["core_vol_median_um3"] = float(np.median(core_vols))

        core_eqsph = eqsph_diams[core_mask]
        result["core_eqsph_diam_mean_um"] = float(np.mean(core_eqsph))
        result["core_eqsph_diam_median_um"] = float(np.median(core_eqsph))

        core_bbox_major = bbox_major[core_mask]
        core_bbox_minor = bbox_minor[core_mask]
        result["core_bbox_major_mean_um"] = float(np.mean(core_bbox_major))
        result["core_bbox_major_median_um"] = float(np.median(core_bbox_major))
        result["core_bbox_minor_mean_um"] = float(np.mean(core_bbox_minor))
        result["core_bbox_minor_median_um"] = float(np.median(core_bbox_minor))

        core_centroids = centroids[core_mask]
        if n_core >= 2:
            core_tree = KDTree(core_centroids)
            core_dists, _ = core_tree.query(core_centroids, k=2)
            core_nn = core_dists[:, 1]
            result["core_nn3d_mean_um"] = float(np.mean(core_nn))
            result["core_nn3d_median_um"] = float(np.median(core_nn))
        else:
            result["core_nn3d_mean_um"] = np.nan
            result["core_nn3d_median_um"] = np.nan

    return result


# ──────────────────────────────────────────────────────────────────────────────
#  Worker for one timepoint (runs in a child process)
# ──────────────────────────────────────────────────────────────────────────────

def _process_one_timepoint(
    tif_path: Path,
    voxel_size: tuple[float, float, float],
) -> tuple[list[dict], dict]:
    """Process a single TIF and return (slice_rows, time_row)."""
    dz, dy, dx = voxel_size
    pixel_area_yx = dy * dx
    tp = _extract_timepoint(tif_path)

    with TiffFile(str(tif_path)) as tif:
        labels_3d = tif.asarray()

    if labels_3d.ndim == 2:
        labels_3d = labels_3d[np.newaxis, ...]

    slice_rows: list[dict] = []
    for z_idx in range(labels_3d.shape[0]):
        s = compute_slice_stats(labels_3d[z_idx], pixel_area_yx)
        s["timepoint"] = tp
        s["z_slice"] = z_idx
        s["z_um"] = z_idx * dz
        slice_rows.append(s)

    # Determine core z-range from per-slice nuclei counts
    n_counts = np.array([s["n_nuclei"] for s in slice_rows])
    peak_n = n_counts.max()
    core_z_range = None
    if peak_n > 0:
        core_slices = np.where(n_counts >= peak_n * CORE_Z_FRAC)[0]
        core_z_range = (int(core_slices.min()), int(core_slices.max()))

    v = compute_volume_stats(labels_3d, voxel_size, core_z_range=core_z_range)
    v["timepoint"] = tp

    return slice_rows, v


# ──────────────────────────────────────────────────────────────────────────────
#  Main pipeline
# ──────────────────────────────────────────────────────────────────────────────

def process_all(
    seg_dir: Path,
    output_dir: Path,
    voxel_override: tuple[float, float, float] | None = None,
    n_workers: int = 1,
) -> None:
    """Process all per-timepoint TIFs and write TSV outputs."""

    tif_files = sorted(
        [f for f in seg_dir.iterdir() if f.suffix == ".tif" and "hyperstack" not in f.stem],
        key=_natural_sort_key,
    )
    if not tif_files:
        print(f"ERROR: No per-timepoint TIF files found in {seg_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(tif_files)} segmented TIFs in {seg_dir}")

    # Read voxel size from first file (or use override)
    if voxel_override:
        voxel_size = voxel_override
        print(f"Using user-supplied voxel size (Z,Y,X): {voxel_size} µm")
    else:
        voxel_size = read_voxel_size(tif_files[0])
        print(f"Auto-detected voxel size (Z,Y,X): {voxel_size} µm")

    dz, dy, dx = voxel_size

    slice_rows: list[dict] = []
    time_rows: list[dict] = []

    n_files = len(tif_files)
    print(f"Processing with {n_workers} worker(s) ...")

    if n_workers == 1:
        # Sequential — simpler progress output
        for i, tif_path in enumerate(tif_files):
            tp = _extract_timepoint(tif_path)
            print(f"  [{i + 1}/{n_files}] t={tp:04d}  {tif_path.name} ...", end="", flush=True)
            s_rows, v = _process_one_timepoint(tif_path, voxel_size)
            slice_rows.extend(s_rows)
            time_rows.append(v)
            print(f"  {v['n_nuclei_3d']} nuclei")
    else:
        # Parallel — submit in sorted order, collect in same order
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = [pool.submit(_process_one_timepoint, p, voxel_size) for p in tif_files]
            for i, (fut, tif_path) in enumerate(zip(futures, tif_files)):
                s_rows, v = fut.result()  # blocks in submission order
                slice_rows.extend(s_rows)
                time_rows.append(v)
                tp = _extract_timepoint(tif_path)
                print(f"  [{i + 1}/{n_files}] t={tp:04d}  {tif_path.name}  {v['n_nuclei_3d']} nuclei")

    # Build DataFrames and sort by timepoint
    df_slice = pd.DataFrame(slice_rows).sort_values(["timepoint", "z_slice"]).reset_index(drop=True)
    df_time = pd.DataFrame(time_rows).sort_values("timepoint").reset_index(drop=True)

    # ── Core-z weighted 2D summary ──────────────────────────────────────────
    # At the embryo cap (first/last z-slices), the imaging plane clips nuclei
    # tangentially, producing under-segmented merged blobs with artifactually
    # large 2D cross-sections.  Rather than hard-trimming, we compute
    # nuclei-count-weighted averages of all per-slice metrics: slices with few
    # nuclei (edges) contribute proportionally less.  We also identify the
    # "core z-range" where the nuclei count exceeds 25 % of the per-timepoint
    # peak — this delineates the reliable interior of the stack.
    # (CORE_Z_FRAC is defined at module level and reused here.)

    core_z_records: list[dict] = []
    for tp, grp in df_slice.groupby("timepoint"):
        n = grp["n_nuclei"].values
        peak_n = n.max()
        if peak_n == 0:
            core_z_records.append({"timepoint": tp})
            continue

        # Core z-range: slices above 25 % of peak nuclei count
        above = grp[grp["n_nuclei"] >= peak_n * CORE_Z_FRAC]
        z_core_lo = int(above["z_slice"].min())
        z_core_hi = int(above["z_slice"].max())

        # Weighted averages over all slices with nuclei
        nz = grp[grp["n_nuclei"] > 0]
        weights = nz["n_nuclei"].values.astype(float)

        def _wmean(col: str) -> float:
            vals = nz[col].values
            mask = ~np.isnan(vals)
            if mask.sum() == 0:
                return np.nan
            return float(np.average(vals[mask], weights=weights[mask]))

        rec: dict = {
            "timepoint": tp,
            "z_core_lo": z_core_lo,
            "z_core_hi": z_core_hi,
            "n_core_slices": z_core_hi - z_core_lo + 1,
            "wt_density_per_um2": _wmean("density_per_um2"),
            "wt_size_mean_um2": _wmean("size_mean_um2"),
            "wt_size_median_um2": _wmean("size_median_um2"),
            "wt_eqcircle_diam_mean_um": _wmean("eqcircle_diam_mean_um"),
            "wt_eqcircle_diam_median_um": _wmean("eqcircle_diam_median_um"),
            "wt_nn_dist_mean_um": _wmean("nn_dist_mean_um"),
            "wt_nn_dist_median_um": _wmean("nn_dist_median_um"),
        }
        core_z_records.append(rec)

    df_core = pd.DataFrame(core_z_records)
    df_time = df_time.merge(df_core, on="timepoint", how="left")

    # Reorder columns
    slice_col_order = [
        "timepoint", "z_slice", "z_um",
        "n_nuclei", "slice_area_um2", "density_per_um2",
        "size_mean_um2", "size_median_um2", "size_std_um2",
        "eqcircle_diam_mean_um", "eqcircle_diam_median_um", "eqcircle_diam_std_um",
        "nn_dist_mean_um", "nn_dist_median_um", "nn_dist_std_um",
    ]
    df_slice = df_slice[[c for c in slice_col_order if c in df_slice.columns]]

    time_col_order = [
        "timepoint",
        "n_nuclei_3d", "volume_um3", "density_per_um3",
        "vol_mean_um3", "vol_median_um3", "vol_std_um3",
        "eqsph_diam_mean_um", "eqsph_diam_median_um", "eqsph_diam_std_um",
        "bbox_major_mean_um", "bbox_major_median_um", "bbox_major_std_um",
        "bbox_minor_mean_um", "bbox_minor_median_um", "bbox_minor_std_um",
        "nn3d_mean_um", "nn3d_median_um", "nn3d_std_um",
        # Core-filtered 3D stats (centroid-z within core range)
        "core_n_nuclei_3d",
        "core_vol_mean_um3", "core_vol_median_um3",
        "core_eqsph_diam_mean_um", "core_eqsph_diam_median_um",
        "core_bbox_major_mean_um", "core_bbox_major_median_um",
        "core_bbox_minor_mean_um", "core_bbox_minor_median_um",
        "core_nn3d_mean_um", "core_nn3d_median_um",
        # 2D summary
        "z_core_lo", "z_core_hi", "n_core_slices",
        "wt_density_per_um2",
        "wt_size_mean_um2", "wt_size_median_um2",
        "wt_eqcircle_diam_mean_um", "wt_eqcircle_diam_median_um",
        "wt_nn_dist_mean_um", "wt_nn_dist_median_um",
    ]
    df_time = df_time[[c for c in time_col_order if c in df_time.columns]]

    # Global summary
    global_stats = {
        "total_timepoints": len(tif_files),
        "total_nuclei_detected": int(df_time["n_nuclei_3d"].sum()),
        "nuclei_per_tp_mean": float(df_time["n_nuclei_3d"].mean()),
        "nuclei_per_tp_std": float(df_time["n_nuclei_3d"].std()),
        "density3d_mean": float(df_time["density_per_um3"].mean()),
        "density3d_std": float(df_time["density_per_um3"].std()),
        "vol_mean_um3": float(df_time["vol_mean_um3"].mean()),
        "eqsph_diam_mean_um": float(df_time["eqsph_diam_mean_um"].mean()),
        "bbox_major_mean_um": float(df_time["bbox_major_mean_um"].mean()),
        "bbox_minor_mean_um": float(df_time["bbox_minor_mean_um"].mean()),
        "nn3d_mean_um": float(df_time["nn3d_mean_um"].mean()),
        "voxel_dz_um": dz,
        "voxel_dy_um": dy,
        "voxel_dx_um": dx,
    }
    df_global = pd.DataFrame([global_stats])

    # Write outputs
    output_dir.mkdir(parents=True, exist_ok=True)

    slice_path = output_dir / "nuclear_stats_per_slice.tsv"
    time_path = output_dir / "nuclear_stats_per_time.tsv"
    global_path = output_dir / "nuclear_stats_global.tsv"

    df_slice.to_csv(slice_path, sep="\t", index=False, float_format="%.6g")
    df_time.to_csv(time_path, sep="\t", index=False, float_format="%.6g")
    df_global.to_csv(global_path, sep="\t", index=False, float_format="%.6g")

    print(f"\nOutputs written to {output_dir}/")
    print(f"  {slice_path.name}  ({len(df_slice)} rows)")
    print(f"  {time_path.name}   ({len(df_time)} rows)")
    print(f"  {global_path.name}     (1 row)")


# ──────────────────────────────────────────────────────────────────────────────
#  CLI
# ──────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute nuclear stats from segmented label TIFs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "segmented_dir",
        type=Path,
        help="Directory with per-timepoint segmented TIFs (t0001_segmented.tif, …)",
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: <segmented_dir>/../nuclear_stats_output)",
    )
    parser.add_argument(
        "--voxel-size",
        nargs=3,
        type=float,
        metavar=("DZ", "DY", "DX"),
        default=None,
        help="Voxel size in µm (Z Y X). Overrides TIFF metadata.",
    )
    parser.add_argument(
        "-j", "--workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: all available CPUs).",
    )
    args = parser.parse_args()

    seg_dir = args.segmented_dir.resolve()
    if not seg_dir.is_dir():
        print(f"ERROR: {seg_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    if args.output_dir:
        out_dir = args.output_dir.resolve()
    else:
        out_dir = seg_dir.parent / "nuclear_stats_output"

    voxel = tuple(args.voxel_size) if args.voxel_size else None
    n_workers = args.workers if args.workers else os.cpu_count() or 1

    process_all(seg_dir, out_dir, voxel_override=voxel, n_workers=n_workers)


if __name__ == "__main__":
    main()
