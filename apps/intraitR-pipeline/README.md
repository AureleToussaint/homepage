# intraitR pipeline — point-and-click interface to the intraitR workflow

A Shiny application that runs the whole [intraitR](https://github.com/FunTraits/intraitR)
pipeline without writing R code:

| Tab | What it does | intraitR / Rfishmorph functions behind it |
|---|---|---|
| 1 · Start | example data, or upload a digitizer workbook (+ optional determinations CSV) | `load_t26_saudrune_landmarks()`, `read_landmarks_xlsx()` |
| 2 · Measure | photographs added as files, as a **whole folder** or as a **ZIP**, then **serves the package's own digitizer** on them and loads what it saved | the Shiny app behind `digitize_landmarks()`, then `consolidate_landmarks()` + `read_landmarks_xlsx()` |
| 3 · Check | per-specimen plot, impute / orientation / geometry corrections, Procrustes outlier screening | `plot_fishmorph_points()`, `impute_landmarks()`, `standardize_orientation()`, `correct_geometry()`, `gpa_fish()`, `detect_outliers()` |
| 4 · Traits | 11 segments (cm) and 9 ratios, Villéger et al. (2010) special cases, summary by species | `fishmorph_segments()`, `fishmorph_ratios()`, `summary_traits()` |
| 5 · Trait space | PCA of the ratios, spider/hull/density display, disparity test | `trait_space()`, `trait_disparity()` |
| 6 · FISHMORPH | projection into the global FISHMORPH morphospace | `project_fishmorph()` |
| 7 · ITV | inter/intraspecific variance partition, CV, accumulation curves | `itv_index()`, `intraspecific_variability()`, `itv_accumulation()` |
| 8 · Shape | GPA + shape-space PCA, shape disparity | `fishmorph_shape_landmarks()`, `gpa_fish()`, `shape_space()`, `intraspecific_variability()` |
| 9 · Repeatability | %ME / repeatability per segment, placement error per landmark (bias sheet) | `measurement_error()`, `digitization_error()` |
| 10 · Export | zip of every table (.csv), figure (.png/.pdf) and the R script equivalent to the session | — |

### Getting the photographs in

Three routes, all appending to the same job (so a big set can go up in several
goes, and stragglers can be added later):

- **Files** — several files at once.
- **Folder** — a whole directory, sub-folders included (`webkitdirectory`:
  Chrome, Edge, Safari; Firefox does not implement it).
- **ZIP** — any zipped folder, flattened on extraction (`junkpaths`, which also
  makes `../` entries harmless).

Non-images are ignored. A name already in the set is **not** renamed: the file
name without its extension is the specimen code, so a silent rename would
invent a specimen — an identical file is skipped, a different one is reported
and left out for you to rename. *Start a new set* empties the session after a
confirmation.

`options(shiny.maxRequestSize = 300 * 1024^2)` raises Shiny's 5 MB default,
which a handful of real photographs would otherwise exceed at once.

### Why the photographs are measured small

The digitizer downsamples every photograph to its *Display* setting (1200 px by
default) the moment it opens it, and drops the original — resolution above that
is uploaded, decoded and held in memory without ever being measured on. Measured
on a 4000 × 3000 frame:

| | 12 Mpx original | 2400 px copy |
|---|---|---|
| opening one specimen | 1.02 s | 0.25 s |
| the decoded array in R | 288 MB | 72 MB |
| one redraw | 0.41 s | 0.15 s |

So photographs are reduced to **2400 px on the long side** (selector in the
sidebar: 2400 / 1600 / original), twice:

- **in the browser**, before the upload, for the *Files* and *Folder* routes —
  `createImageBitmap(..., imageOrientation: "from-image")` (so the EXIF rotation
  is applied, which a bare canvas drops) then a canvas and `toBlob`, with the
  resized `File` objects put back on the input through a `DataTransfer`. This is
  what cuts the upload: a 5 MB frame leaves as ~200 kB. A browser that refuses
  any of it simply uploads the original.
- **on the server**, for anything that still arrives larger — a ZIP (the archive
  is opened server-side, so the browser cannot shrink it), an old browser, a
  failed canvas. `magick` does it with EXIF rotation and proper resampling; with
  magick absent, a base-R subsample takes over, which is exactly what the
  digitizer itself does to display a photograph.

**This changes no measured trait.** The nine FISHMORPH ratios are quotients of
segments of the same image, and the scale bar is digitized on that same image,
so a uniform resize cancels out. What it costs is placement precision, and
2400 px leaves ~0.04 mm/px on a 10 cm fish, two orders of magnitude finer than
the ~0.3 % of standard length that the package's own T-26 repeat trial measures
as digitization bias. A ZIP of already-small photographs uploads fastest of all.

### What the pipeline tab does while you measure

Both tabs share one R process on shinyapps, so anything the pipeline tab does
while the digitizer is open is time the digitizer does not get. The journal is
therefore re-read only when its files have actually changed (a signature built
from the directory listing), only every 10 s, and only while the *Measure* tab
is the one on screen.

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

## Live at https://globaltrait.shinyapps.io/intraitR-pipeline/

That is the URL the *Open the app* button of the site's R packages page points to
(`Rpackages.Rmd`, card "Pipeline app", and the students' block above it).
Redeploying with the same `appName` replaces it in place; the link does not change.

## Deploying (shinyapps.io, free tier is enough to start)

1. Create an account on https://www.shinyapps.io and, in RStudio, paste the
   token from *Account → Tokens* (`rsconnect::setAccountInfo(...)`).
2. Install the packages above **from their sources** (`remotes::install_github()`
   for Rfishmorph and intraitR): rsconnect records the GitHub origin and
   reinstalls them on the server.
3. `source("deploy.R")` — or `rsconnect::deployApp("intraitR-pipeline", appName = "intraitR-pipeline")`.
4. The app is then at `https://<account>.shinyapps.io/intraitR-pipeline/` —
   here `globaltrait`. A different account or `appName` means a different URL, so
   update the two links in `Rpackages.Rmd` (and the rendered `Rpackages.html`).

Notes: the first launch loads the FISHMORPH reference (~9 000 species) when the
FISHMORPH tab is used; the free tier (1 GB RAM, 25 active hours / month) copes
with a class of students, but set *Instance size* to "Large" in the shinyapps
dashboard if projections fail. `itv_accumulation()` with many permutations is
the slowest step; the default (49) keeps it under a minute.

## Files

- `app.R` — the application (single file).
- `deploy.R` — deployment helper for shinyapps.io.
- `README.md` — this file.
