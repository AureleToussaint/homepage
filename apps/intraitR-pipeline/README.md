# intraitR pipeline — point-and-click interface to the intraitR workflow

A Shiny application that runs the whole [intraitR](https://github.com/FunTraits/intraitR)
pipeline without writing R code:

| Tab | What it does | intraitR / Rfishmorph functions behind it |
|---|---|---|
| 1 · Start | example data, or upload a digitizer workbook (+ optional determinations CSV) | `load_t26_saudrune_landmarks()`, `read_landmarks_xlsx()` |
| 2 · Measure | **serves the package's own digitizer** on the uploaded photographs, then loads what it saved | the Shiny app behind `digitize_landmarks()`, then `consolidate_landmarks()` + `read_landmarks_xlsx()` |
| 3 · Check | per-specimen plot, impute / orientation / geometry corrections, Procrustes outlier screening | `plot_fishmorph_points()`, `impute_landmarks()`, `standardize_orientation()`, `correct_geometry()`, `gpa_fish()`, `detect_outliers()` |
| 4 · Traits | 11 segments (cm) and 9 ratios, Villéger et al. (2010) special cases, summary by species | `fishmorph_segments()`, `fishmorph_ratios()`, `summary_traits()` |
| 5 · Trait space | PCA of the ratios, spider/hull/density display, disparity test | `trait_space()`, `trait_disparity()` |
| 6 · FISHMORPH | projection into the global FISHMORPH morphospace | `project_fishmorph()` |
| 7 · ITV | inter/intraspecific variance partition, CV, accumulation curves | `itv_index()`, `intraspecific_variability()`, `itv_accumulation()` |
| 8 · Shape | GPA + shape-space PCA, shape disparity | `fishmorph_shape_landmarks()`, `gpa_fish()`, `shape_space()`, `intraspecific_variability()` |
| 9 · Repeatability | %ME / repeatability per segment, placement error per landmark (bias sheet) | `measurement_error()`, `digitization_error()` |
| 10 · Export | zip of every table (.csv), figure (.png/.pdf) and the R script equivalent to the session | — |

## The *Measure* tab is the package's digitizer, not a copy of it

`digitize_landmarks()` hands its session configuration to its Shiny application
in one option (`intraitR.digitizer`), read when the app file is evaluated. The
pipeline app does the same: it writes the uploaded photographs to a job folder,
sets that option, sources
`system.file("shiny/landmarking_app", package = "intraitR")` into a fresh
environment, and serves its `ui`/`server` under `?page=digitizer&job=<token>` —
so the digitizer opens in a second browser tab as its own full page, with its
queues (New / Correct / Repeats), its landmark buttons and statuses, its zoom
and pan, its protocol checks (eye vertical, extreme points, coincidence rules),
its plates and its repeat mode. Nothing is reimplemented, so the app cannot
drift from the package.

What it writes is therefore what the package writes: one workbook plus the
append-only journal. *Load the measurements into the pipeline* rebuilds the
table from the journal with `consolidate_landmarks()` and reads it back with
`read_landmarks_xlsx()` — the exact route a campaign takes — and the zip from
*Download workbook + journal* is what a student sends back for merging.

The only difference from a laptop: landmark **prediction** needs the ml-morph
model and a Python/dlib environment, which a server does not have, so the
*Predict 19 landmarks* button reports that no model is available and every
point is placed by hand. To offer prediction online, install the model
directory on the server and set `INTRAITR_MLMORPH_DIR` before the app starts.

## Run locally

```r
install.packages(c("shiny", "bslib", "geomorph", "readxl", "writexl", "jpeg", "png", "zip", "remotes"))
remotes::install_github("FunTraits/Rfishmorph")
remotes::install_github("FunTraits/intraitR")
shiny::runApp("intraitR-pipeline")
```

## Deploy on shinyapps.io (free tier is enough to start)

1. Create an account on https://www.shinyapps.io and, in RStudio, paste the
   token from *Account → Tokens* (`rsconnect::setAccountInfo(...)`).
2. Install the packages above **from their sources** (`remotes::install_github()`
   for Rfishmorph and intraitR): rsconnect records the GitHub origin and
   reinstalls them on the server.
3. `source("deploy.R")` — or `rsconnect::deployApp("intraitR-pipeline", appName = "intraitR-pipeline")`.
4. The app is then at `https://<account>.shinyapps.io/intraitR-pipeline/`;
   put that URL in `Rpackages.Rmd` (card "Pipeline app").

Notes: the first launch loads the FISHMORPH reference (~9 000 species) when the
FISHMORPH tab is used; the free tier (1 GB RAM, 25 active hours / month) copes
with a class of students, but set *Instance size* to "Large" in the shinyapps
dashboard if projections fail. `itv_accumulation()` with many permutations is
the slowest step; the default (49) keeps it under a minute.

## Files

- `app.R` — the application (single file).
- `deploy.R` — deployment helper for shinyapps.io.
- `README.md` — this file.
