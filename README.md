# SPIM Medaka / Zebrafish Trajectory Analysis

Comparative 4D lightsheet whole-embryo nuclei tracking analysis — **medaka** vs **zebrafish** gastrulation, with a focus on cell ingression at the margin.

## Pipeline overview

| Step | Script | Language | Purpose |
|------|--------|----------|---------|
| 1 | `embryo_viewer.py` | Python (napari) | Interactive 4D viewer: orientation, sphere fit, margin selection, export |
| 1′ | `orientation_interactive.R` | R (Shiny) | Alternative Step 1 (lighter weight, no GPU needed) |
| 2 | `trackmate_filter_and_validate.R` | R | Parameter-sweep QC filtering against manual ground truth |
| 3 | `trackmate_analysis.R` | R | Full biological analysis (speed, displacement, ingression, radial velocity, …) |

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

# Creates .venv/ and installs all pinned dependencies (napari, scipy, pandas, …)
uv sync
```

> `uv sync` reads `pyproject.toml` + `uv.lock` and reproduces the exact same
> environment on any machine. No conda, no Docker needed.

### 2. Set up R environment

```r
# In R, from the project directory:
renv::restore()
```

### 3. Place your data

Put TrackMate CSV exports in the project root:

```
medaka_25082025_combined_spots.csv
medaka_25082025_combined_tracks.csv
spots_zebrafish.csv
tracks_zebrafish.csv
```

*(CSVs are git-ignored because they are too large.)*

## Running the napari 4D viewer (Step 1)

```bash
# Full dataset
uv run python embryo_viewer.py medaka_25082025_combined_spots.csv

# With tracks overlay (progressive formation as time advances)
uv run python embryo_viewer.py medaka_25082025_combined_spots.csv \
    -t medaka_25082025_combined_tracks.csv

# Subsample for faster exploration on a laptop
uv run python embryo_viewer.py medaka_25082025_combined_spots.csv --max-spots 500000
```

**Inside the viewer:**

1. Use the **time slider** (bottom) to scrub through frames — tracks form progressively
2. **Pick Animal Pole** → click a nucleus → **Pick Dorsal** → click a nucleus
3. **Orient Embryo** → standardises axes (AP → +Y, Dorsal → +X)
4. **Fit Sphere** → cap-aware sphere fit for SPIM thin-slab data
5. **Apply Margin** → flag spots within a latitude band (adjustable sliders)
6. **Colour buttons** → depth / latitude / frame / margin overlay
7. **Export** → writes `analysis_output/oriented_spots.csv` + `sphere_params.csv`

### Running the R Shiny alternative (Step 1′)

```r
# In RStudio or R console
shiny::runApp("orientation_interactive.R")
```

### Running Steps 2 & 3

```r
source("trackmate_filter_and_validate.R")
source("trackmate_analysis.R")
```

## Project structure

```
├── embryo_viewer.py               # napari 4D viewer (Python)
├── orientation_interactive.R      # Shiny orientation app (R, alternative)
├── trackmate_filter_and_validate.R
├── trackmate_analysis.R
├── pyproject.toml                 # Python dependencies (source of truth)
├── uv.lock                        # Pinned Python lockfile (commit this!)
├── requirements.txt               # Pip fallback
├── renv.lock                      # Pinned R lockfile
├── renv/activate.R
├── .gitignore
└── archive/                       # Previous script versions
```

## Key files to commit

```
# Code
embryo_viewer.py
orientation_interactive.R
trackmate_filter_and_validate.R
trackmate_analysis.R

# Environment lockfiles (reproducibility)
pyproject.toml
uv.lock
requirements.txt
renv.lock
renv/activate.R
renv/settings.json

# Project config
.gitignore
spim_medaka_zebrafish_trajectory.Rproj
README.md
```

## Reproducibility notes

- **Python**: `uv.lock` pins every transitive dependency to an exact version + hash.
  Anyone running `uv sync` gets the identical environment.
- **R**: `renv.lock` does the same for R packages. `renv::restore()` reproduces it.
- **Data**: CSVs are git-ignored. Store them on a shared drive / Zenodo / figshare
  and document the path in your lab notebook.
