# SPIM Medaka / Zebrafish Trajectory Analysis

Comparative 4D lightsheet whole-embryo nuclei tracking analysis — **medaka** vs **zebrafish** gastrulation, with a focus on cell ingression at the margin.

## Pipeline overview

| Step | Script | Language | Purpose |
|------|--------|----------|---------|
| 1 | `embryo_viewer.py` | Python (napari) | Interactive 4D viewer: orientation, sphere fit, margin selection, export |
| 1′ | `ultrack_viewer.py` | Python (napari) | Variant of viewer for ultrack-format data |
| 2 | `nuclear_stats.py` | Python | Per-timepoint nuclear density, size, and internuclear distance from segmented labels |
| 3 | `gastrulation_dynamics_medaka.R` | R | Medaka-specific dynamics analysis (10 figures) |
| 3′ | `gastrulation_dynamics_zebrafish.R` | R | Zebrafish-specific dynamics analysis (10 figures) |
| 4 | `gastrulation_dynamics_comparison.R` | R | Cross-species comparison (7 figures) |
| S1 | `medaka_density_analysis.R` | R | Supplementary: density & flow figures (medaka) |
| S2 | `zebrafish_density_analysis.R` | R | Supplementary: density & flow figures (zebrafish) |

## Quick start

### Prerequisites

| Tool | Install |
|------|---------|
| **uv** (Python env manager) | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |
| **R ≥ 4.5** + renv | Already configured in `renv.lock` |

### 1. Clone & set up Python environment

```bash
git clone <repo-url>
cd spim_medaka_zebrafish_trajectory
uv sync
```

### 2. Set up R environment

```r
renv::restore()
```

### 3. Place your data

Put TrackMate CSV exports in the project root:

```
medaka_25082025_combined_spots.csv
medaka_25082025_combined_tracks.csv
```

*(CSVs are git-ignored because they are too large.)*

## Running the napari 4D viewer (Step 1)

```bash
uv run python embryo_viewer.py medaka_25082025_combined_spots.csv \
    -t medaka_25082025_combined_tracks.csv
```

**Inside the viewer:**

1. Use the **time slider** to scrub through frames — tracks form progressively
2. **Pick Animal Pole** → click a nucleus → **Pick Dorsal** → click a nucleus
3. **Orient Embryo** → standardises axes (AP → +Y, Dorsal → +X)
4. **Fit Sphere** → cap-aware sphere fit for SPIM thin-slab data
5. **Apply Margin** → flag spots within a latitude band (adjustable sliders)
6. **Colour buttons** → depth / latitude / frame / margin overlay
7. **Export** → writes oriented CSVs + `sphere_params.csv` to `oriented_*_ultrack/`

## Running the R analysis (Steps 3–4)

```r
source("gastrulation_dynamics_medaka.R")     # → analysis_output_medaka/
source("gastrulation_dynamics_zebrafish.R")   # → analysis_output_zebrafish/
source("gastrulation_dynamics_comparison.R")  # → analysis_output_comparison/
```

## Project structure

```
├── embryo_viewer.py                     # napari 4D viewer (Step 1)
├── ultrack_viewer.py                    # ultrack-format variant
├── nuclear_stats.py                     # nuclear segmentation stats (Step 2)
├── gastrulation_dynamics_medaka.R       # medaka analysis (Step 3)
├── gastrulation_dynamics_zebrafish.R    # zebrafish analysis (Step 3′)
├── gastrulation_dynamics_comparison.R   # cross-species comparison (Step 4)
├── medaka_density_analysis.R            # supplementary density figures
├── zebrafish_density_analysis.R         # supplementary density figures
├── gastrulation_dynamics_methods.tex    # methods write-up (LaTeX)
│
├── oriented_medaka_ultrack/             # oriented input data (medaka)
├── oriented_zebrafish_ultrack/          # oriented input data (zebrafish)
├── nuclei_stats_medaka/                 # nuclear stats TSVs (medaka)
├── nuclei_stats_zebrafish/              # nuclear stats TSVs (zebrafish)
├── medaka_25082025_combined_*.csv       # raw TrackMate exports (git-ignored)
│
├── analysis_output_medaka/              # medaka output (10 PDFs + CSVs)
├── analysis_output_zebrafish/           # zebrafish output (10 PDFs + CSVs)
├── analysis_output_comparison/          # comparison output (7 PDFs + CSV)
├── analysis_output_medaka_25082025/     # density analysis input (medaka)
├── analysis_output_zebrafish_05112025/  # density analysis input (zebrafish)
│
├── pyproject.toml                       # Python dependencies
├── uv.lock                             # pinned Python lockfile
├── requirements.txt                     # pip fallback
├── renv.lock                            # pinned R lockfile
├── renv/                                # R environment
├── .gitignore
├── README.md
└── archive/                             # superseded scripts & data
```

## Reproducibility notes

- **Python**: `uv sync` reproduces the exact environment from `pyproject.toml` + `uv.lock`.
- **R**: `renv::restore()` reproduces the R environment from `renv.lock`.
- **Data**: CSVs are git-ignored. Store them on a shared drive / Zenodo / figshare.
