# =============================================================================
# intraitR pipeline -- a point-and-click interface to the intraitR workflow
#
#   photographs -> landmarks -> quality control -> FISHMORPH traits ->
#   trait space / FISHMORPH morphospace / intraspecific variability /
#   shape space / repeatability -> figures, tables and an R script
#
# Every computation is done by an exported intraitR (or Rfishmorph) function;
# the app only collects the arguments, calls the function and shows the result.
# The measuring step is not a re-implementation either: it SERVES the package's
# own digitizer (the Shiny application behind intraitR::digitize_landmarks()),
# configured per job and opened in a second browser tab, so that a student and
# a user of the package see exactly the same interface and produce exactly the
# same workbook and append-only journal.
# The R code equivalent to what was clicked is accumulated and can be exported
# (Export tab), so a user who wants to move to R starts from a working script.
#
# Run locally:      shiny::runApp("intraitR-pipeline")
# Deploy:           rsconnect::deployApp("intraitR-pipeline")   (see README.md)
#
# Author: Aurele Toussaint (CNRS, CRBE) -- https://github.com/FunTraits
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(intraitR)   # Rfishmorph is imported by intraitR; it is not attached on purpose:
                      # both packages export fishmorph_segments()/fishmorph_ratios() and the
                      # intraitR versions (scale_action, na_action, ...) are the ones used here.
})

APP_VERSION <- "0.2.0"
RATIO_COLS  <- c("BEl", "VEp", "REs", "OGp", "RMl", "BLs", "PFv", "PFs", "CPt")
SEG_COLS    <- c("Bl", "Bd", "Hd", "Eh", "Mo", "PFi", "PFl", "Ed", "Jl", "CPd", "CFd")
N_LM_READ   <- 23L   # what an analysis reads back from a digitizer workbook

has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

# ----------------------------------------------------------------------------- the digitizer
# The measuring step is the package's own application, not a copy of it:
# `system.file("shiny/landmarking_app", package = "intraitR")` is the very file
# that intraitR::digitize_landmarks() runs. That launcher hands its session
# configuration over in ONE option (`intraitR.digitizer`), read when the app
# file is evaluated -- so a job is opened by setting the option and sourcing the
# file into a fresh environment, whose `ui` and `server` are then served under
# ?page=digitizer&job=<token>. Each job keeps its own photographs, workbook and
# journal; the environment holds no mutable state (everything a session touches
# lives inside server()), so several browser tabs may share one job safely.
DIGITIZER_DIR <- function() system.file("shiny", "landmarking_app", package = "intraitR")

JOBS <- new.env(parent = emptyenv())

job_new <- function(token, cfg) assign(token, list(cfg = cfg, env = NULL), envir = JOBS)

job_get <- function(token) if (!is.null(token) && exists(token, envir = JOBS, inherits = FALSE))
  get(token, envir = JOBS, inherits = FALSE) else NULL

# Build (once per job) the environment the digitizer app lives in.
job_env <- function(token) {
  job <- job_get(token)
  if (is.null(job)) return(NULL)
  if (!is.null(job$env)) return(job$env)
  dir <- DIGITIZER_DIR()
  if (!nzchar(dir) || !file.exists(file.path(dir, "app.R"))) return(NULL)
  old <- options(intraitR.digitizer = job$cfg)
  on.exit(options(old), add = TRUE)
  e <- new.env(parent = globalenv())
  # Read as UTF-8 and parse as UTF-8, the way shiny::runApp() reads an app file:
  # the digitizer's labels are full of em dashes, and source(encoding = "UTF-8")
  # re-encodes to the native locale, which fails outright under a C locale.
  ok <- tryCatch({
    lines <- readLines(file.path(dir, "app.R"), warn = FALSE, encoding = "UTF-8")
    for (ex in parse(text = lines, keep.source = FALSE, encoding = "UTF-8")) eval(ex, envir = e)
    TRUE
  }, error = function(err) { warning("digitizer: ", conditionMessage(err)); FALSE })
  if (!ok || is.null(e$ui) || !is.function(e$server)) return(NULL)
  job$env <- e; assign(token, job, envir = JOBS)
  e
}

job_token_of <- function(search) {
  q <- shiny::parseQueryString(search %||% "")
  if (identical(q$page, "digitizer")) q$job else NULL
}

PHOTO_EXT <- c("jpg", "jpeg", "JPG", "JPEG", "png", "PNG", "tif", "tiff", "bmp", "gif")

# The photograph a specimen code came from (plates add an _i<k> suffix).
photo_for_code <- function(dir, code) {
  if (is.null(dir) || !nzchar(dir) || !dir.exists(dir)) return(NULL)
  base <- sub("_i[0-9]+$", "", code)
  cand <- file.path(dir, paste0(base, ".", PHOTO_EXT))
  cand <- cand[file.exists(cand)]
  if (length(cand)) cand[1] else NULL
}


# ----------------------------------------------------------------------------- helpers
`%||%` <- function(a, b) if (is.null(a)) b else a

notify_error <- function(e, where = "") {
  showNotification(paste0(if (nzchar(where)) paste0(where, ": ") else "", conditionMessage(e)),
                   type = "error", duration = 12)
}

# Run an expression, show a progress bar, return NULL (and a notification) on error
safely <- function(expr, where = "", message = "Computing...") {
  withProgress(message = message, value = 0.3, {
    tryCatch(expr, error = function(e) { notify_error(e, where); NULL })
  })
}

step_card <- function(title, ..., icon = NULL) {
  card(card_header(class = "bg-light fw-semibold", if (!is.null(icon)) icon(icon), " ", title), ...)
}

help_text <- function(...) tags$p(class = "text-muted small mb-2", ...)

ratio_table <- function(df) {
  df[] <- lapply(df, function(x) {
    if (!is.numeric(x)) return(x)
    if (all(is.na(x) | x == round(x))) return(as.integer(round(x)))   # counts stay integers
    round(x, 4)
  })
  df
}

detect_layout <- function(cols) {
  if (any(grepl("^[0-9]+_X$", cols))) {
    n <- sum(grepl("^[0-9]+_X$", cols))
    list(x = "{i}_X", y = "{i}_Y", n = min(n, N_LM_READ), kind = "digitizer")
  } else if (any(grepl("^X_[0-9]+$", cols))) {
    list(x = "X_{i}", y = "Y_{i}", n = sum(grepl("^X_[0-9]+$", cols)), kind = "classic")
  } else NULL
}

# the species vector of a landmark object (NULL if none)
species_of <- function(lm) {
  if (is.null(lm) || is.null(lm$metadata)) return(NULL)
  s <- lm$metadata$species
  if (is.null(s) || all(is.na(s))) return(NULL)
  s
}

# ----------------------------------------------------------------------------- UI
theme <- bs_theme(version = 5, bootswatch = "cosmo", primary = "#27ae60",
                  "navbar-bg" = "#2c3e50")

sidebar_w <- 330

pipeline_ui <- page_navbar(
  title = tags$span(tags$b("intraitR"), " pipeline ", tags$small(class = "opacity-75", paste0("v", APP_VERSION))),
  theme = theme, id = "nav", fillable = FALSE,
  header = tags$head(tags$style(HTML("
    .lm-badge { display:inline-block; width:30px; height:30px; border-radius:50%; background:#27ae60;
                color:#fff; text-align:center; line-height:30px; font-weight:700; margin-right:8px; }
    .status-ok { color:#27ae60; font-weight:600; } .status-no { color:#999; }
    pre.rcode { background:#f7f9fb; border:1px solid #e4e9f0; border-radius:8px; padding:10px; font-size:.82em; }
    .nav-link { font-weight:500; }
    .card { margin-bottom: 14px; }
  "))),

  # ---------------------------------------------------------------- 1. Start
  nav_panel("1 · Start", icon = icon("play"),
    layout_sidebar(width = sidebar_w,
      sidebar = sidebar(
        tags$h5("Where do your data come from?"),
        help_text("Pick one of the three options. You can always come back here and start again."),
        tags$hr(),
        tags$b("A. Try with the example data"),
        help_text("The T-26 La Saudrune electrofishing survey shipped with intraitR (7 species, real field photographs already digitized)."),
        actionButton("load_example", "Load the example data set", class = "btn-primary w-100"),
        tags$hr(),
        tags$b("B. Upload a workbook of landmarks"),
        help_text("The .xlsx written by the digitizer (sheet 'measurements', columns 1_X, 1_Y, ...) or any wide sheet with X_1, Y_1, ... columns."),
        fileInput("upload_xlsx", NULL, accept = c(".xlsx", ".xls"), buttonLabel = "Browse...", placeholder = "landmarks.xlsx"),
        uiOutput("sheet_ui"),
        textInput("id_col", "Specimen column", value = "specimen"),
        fileInput("upload_det", "Optional: determinations table (CSV with 'uid' and 'species')", accept = c(".csv", ".txt"),
                  buttonLabel = "Browse...", placeholder = "determinations.csv"),
        actionButton("load_xlsx", "Load the workbook", class = "btn-primary w-100"),
        tags$hr(),
        tags$b("C. Measure photographs yourself"),
        help_text("No landmarks yet? The 'Measure' tab opens intraitR's own digitizer on your photographs."),
        actionButton("goto_measure", "Go to Measure", class = "btn-outline-primary w-100")
      ),
      step_card("Current data set", icon = "database",
        uiOutput("status_box")
      ),
      step_card("How this app works", icon = "circle-info",
        tags$ol(
          tags$li(tags$b("Start / Measure"), " -- get landmark coordinates in: the example, a workbook, or your own photographs digitized with intraitR's own digitizer."),
          tags$li(tags$b("Check"), " -- look at each specimen, impute a missing point, fix orientation and geometry, screen outliers."),
          tags$li(tags$b("Traits"), " -- the 11 FISHMORPH segments (cm) and the 9 dimensionless ratios."),
          tags$li(tags$b("Trait space / FISHMORPH / ITV / Shape"), " -- the figures and tables, each with a few options."),
          tags$li(tags$b("Repeatability"), " -- how reproducible your own clicks are (needs repeated digitizations)."),
          tags$li(tags$b("Export"), " -- every table, every figure and the R script that reproduces them.")
        ),
        help_text("All computations are intraitR / Rfishmorph functions (Toussaint 2026, ",
                  tags$a(href = "https://github.com/FunTraits/intraitR", target = "_blank", "github.com/FunTraits/intraitR"),
                  "); the FISHMORPH protocol follows Brosse et al. (2021).")
      )
    )
  ),

  # ---------------------------------------------------------------- 2. Measure
  nav_panel("2 · Measure", icon = icon("crosshairs"),
    layout_sidebar(width = sidebar_w,
      sidebar = sidebar(
        tags$h5("Digitize your photographs"),
        help_text("This step opens the digitizer of the intraitR package itself — the same application as ",
                  tags$code("digitize_landmarks()"), ", with its queues, its landmark buttons, its zoom and its checks. ",
                  "It writes one workbook and an append-only journal, exactly as it does on a laptop."),
        fileInput("photos", "Photographs (.jpg / .png)", multiple = TRUE,
                  accept = c(".jpg", ".jpeg", ".png", ".tif", ".tiff", ".bmp", ".gif"),
                  buttonLabel = "Browse...", placeholder = "no file selected"),
        help_text("One file per photograph; the file name without its extension is the specimen code. A photograph may hold several fish (a plate): say so in the digitizer."),
        textInput("operator", "Your operator code (initials)", value = ""),
        numericInput("ruler_mm", "Scale bar 20-21: real length (mm)", value = 10, min = 0.1, step = 1),
        numericInput("ind_per_photo", "Fish per photograph", value = 1, min = 1, step = 1),
        selectInput("digit_mode", "Starting queue",
                    choices = c("New photographs" = "new", "Review what is digitized" = "correct",
                                "Repeats (measurement error)" = "repeat")),
        tags$hr(),
        uiOutput("digit_open_ui"),
        tags$hr(),
        actionButton("digit_load", "Load the measurements into the pipeline", class = "btn-primary w-100", icon = icon("arrow-down")),
        help_text("Reads the journal of this session (the source of truth) and brings every saved specimen in. Click it again after measuring more."),
        fileInput("digit_det", "Optional: determinations table (CSV with 'uid' and 'species')",
                  accept = c(".csv", ".txt"), buttonLabel = "Browse...", placeholder = "determinations.csv"),
        downloadButton("digit_download", "Download workbook + journal (.zip)", class = "btn-outline-primary w-100")
      ),
      step_card("How it works", icon = "circle-info",
        tags$ol(
          tags$li("Upload your photographs on the left, give your initials and the real length of the scale bar."),
          tags$li(tags$b("Open the digitizer"), " — it opens in a second browser tab. Place the landmarks there and press ",
                  tags$b("Save & next"), " for each specimen (the queue, the buttons 1–25, the zoom and the checks are the package's own)."),
          tags$li("Come back to this tab and click ", tags$b("Load the measurements into the pipeline"), " — then carry on with tabs 3 to 9."),
          tags$li("Keep ", tags$b("Download workbook + journal"), ": that zip is what you send me, and what gets merged into the campaign.")
        ),
        help_text("Prediction of the landmarks by the ml-morph model is only available when the app runs on a machine where that model and its Python environment are installed; online, every point is placed by hand. Everything else is identical.")
      ),
      step_card("This digitizing session", icon = "folder-open",
        verbatimTextOutput("digit_status"),
        div(style = "overflow-x:auto;", tableOutput("digit_progress"))
      )
    )
  ),

  # ---------------------------------------------------------------- 3. Check
  nav_panel("3 · Check", icon = icon("magnifying-glass"),
    layout_sidebar(width = sidebar_w,
      sidebar = sidebar(
        tags$h5("Quality control"),
        help_text("Look at the specimens one by one, then apply the corrections below. Every correction is recorded in the object (audit trail) and highlighted on the plot: red = imputed point, blue = corrected point."),
        uiOutput("qc_specimen_ui"),
        checkboxInput("qc_photo", "Show the photograph behind the points (when available)", value = TRUE),
        tags$hr(),
        checkboxInput("qc_impute", "Estimate missing anatomical landmarks (impute_landmarks)", value = TRUE),
        checkboxInput("qc_orient", "Enforce head-left / belly-down orientation (standardize_orientation)", value = TRUE),
        checkboxInput("qc_geom", "Standardize geometry and snap to the FISHMORPH conventions (correct_geometry)", value = TRUE),
        actionButton("qc_apply", "Apply corrections", class = "btn-primary w-100", icon = icon("wand-magic-sparkles")),
        tags$hr(),
        numericInput("qc_thr", "Outlier threshold (median + k x MAD)", value = 3, min = 1, step = 0.5),
        actionButton("qc_outliers", "Screen shape outliers (gpa_fish + detect_outliers)", class = "btn-outline-primary w-100")
      ),
      layout_columns(col_widths = c(7, 5),
        step_card("Specimen", icon = "fish",
          plotOutput("qc_plot", height = "520px")
        ),
        step_card("Corrections and outliers", icon = "clipboard-check",
          verbatimTextOutput("qc_summary"),
          div(style = "overflow-x:auto; font-size:.85em;", tableOutput("qc_outlier_table"))
        )
      )
    )
  ),

  # ---------------------------------------------------------------- 4. Traits
  nav_panel("4 · Traits", icon = icon("ruler"),
    layout_sidebar(width = sidebar_w,
      sidebar = sidebar(
        tags$h5("FISHMORPH traits"),
        help_text("11 linear segments converted to centimetres with the scale bar (points 20-21), then the 9 dimensionless ratios of Villeger et al. (2010) / Brosse et al. (2021)."),
        numericInput("scale_cm", "Real length of the scale bar (cm)", value = 1, min = 0.01, step = 0.1),
        radioButtons("scale_action", "Specimen without a usable scale bar",
                     choices = c("Segments set to NA (ratios unaffected)" = "na", "Segments kept in pixels" = "pixels")),
        radioButtons("na_action_tr", "Missing values in the ratios",
                     choices = c("Keep NA" = "keep", "Drop the specimen" = "omit", "Impute species mean" = "impute_group_mean")),
        tags$b("Special morphologies (Villeger et al. 2010)"),
        checkboxInput("no_caudal", "No visible caudal fin (CPt = 1)", FALSE),
        checkboxInput("ventral_mouth", "Ventral mouth, algae browser (OGp = RMl = 0)", FALSE),
        checkboxInput("no_pectoral", "No pectoral fin (PFv = 0)", FALSE),
        actionButton("traits_run", "Compute the traits", class = "btn-primary w-100", icon = icon("calculator")),
        tags$hr(),
        downloadButton("dl_segments", "Segments (.csv)", class = "w-100 mb-1"),
        downloadButton("dl_ratios", "Ratios (.csv)", class = "w-100 mb-1"),
        downloadButton("dl_summary", "Summary by species (.csv)", class = "w-100")
      ),
      step_card("Ratios (one row per specimen)", icon = "table",
        verbatimTextOutput("traits_msg"),
        div(style = "overflow-x:auto;", tableOutput("ratios_table"))
      ),
      step_card("Summary by species (mean, sd, min, max)", icon = "chart-simple",
        div(style = "overflow-x:auto;", tableOutput("summary_table"))
      )
    )
  ),

  # ---------------------------------------------------------------- 5. Trait space
  nav_panel("5 · Trait space", icon = icon("chart-area"),
    layout_sidebar(width = sidebar_w,
      sidebar = sidebar(
        tags$h5("Functional trait space of your sample"),
        help_text("PCA of the 9 log10(x+1)-transformed, standardised ratios (trait_space). Each species is a group."),
        selectInput("ts_style", "Display", choices = c("Spider + 95% ellipse" = "spider", "Convex hull" = "hull", "Density contour" = "density", "Points only" = "none")),
        selectInput("ts_axes", "Axes", choices = c("PC1 x PC2" = "1,2", "PC1 x PC3" = "1,3", "PC2 x PC3" = "2,3")),
        checkboxInput("ts_outliers", "Remove within-species outliers before the ordination", FALSE),
        checkboxInput("ts_italic", "Italic species names", TRUE),
        actionButton("ts_run", "Build the trait space", class = "btn-primary w-100", icon = icon("chart-area")),
        tags$hr(),
        numericInput("td_iter", "Permutations for the disparity test", value = 199, min = 99, step = 100),
        actionButton("td_run", "Test differences in trait dispersion (trait_disparity)", class = "btn-outline-primary w-100"),
        tags$hr(),
        downloadButton("dl_ts_png", "Figure (.png)", class = "w-100 mb-1"),
        downloadButton("dl_ts_pdf", "Figure (.pdf)", class = "w-100")
      ),
      step_card("Trait space", icon = "chart-area",
        plotOutput("ts_plot", height = "560px"),
        verbatimTextOutput("ts_print")
      ),
      step_card("Trait disparity (variance of traits per species, permutation test)", icon = "scale-balanced",
        verbatimTextOutput("td_print")
      )
    )
  ),

  # ---------------------------------------------------------------- 6. FISHMORPH
  nav_panel("6 · FISHMORPH", icon = icon("globe"),
    layout_sidebar(width = sidebar_w,
      sidebar = sidebar(
        tags$h5("Your specimens in the global morphospace"),
        help_text("Projection of your individuals into the PCA of the FISHMORPH database (~9 000 freshwater fish species, Brosse et al. 2021) with project_fishmorph. The grey background is the density of the database."),
        selectInput("pf_style", "Display", choices = c("Convex hull per species" = "hull", "Spider" = "spider", "Density" = "density", "Points" = "points")),
        checkboxInput("pf_itvref", "Show the reference position of each species (itv_reference)", TRUE),
        checkboxInput("pf_arrows", "Show the trait loadings (arrows)", FALSE),
        actionButton("pf_run", "Project into FISHMORPH", class = "btn-primary w-100", icon = icon("globe")),
        tags$hr(),
        downloadButton("dl_pf_png", "Figure (.png)", class = "w-100 mb-1"),
        downloadButton("dl_pf_pdf", "Figure (.pdf)", class = "w-100")
      ),
      step_card("FISHMORPH morphospace", icon = "globe",
        plotOutput("pf_plot", height = "600px"),
        verbatimTextOutput("pf_print")
      )
    )
  ),

  # ---------------------------------------------------------------- 7. ITV
  nav_panel("7 · ITV", icon = icon("arrows-left-right"),
    layout_sidebar(width = sidebar_w,
      sidebar = sidebar(
        tags$h5("Intraspecific variability"),
        help_text("How much of the trait variance lies within species (ITV) versus between species (itv_index), the coefficient of variation per species and trait (intraspecific_variability), and how the within-species variance accumulates with the number of individuals (itv_accumulation)."),
        actionButton("itv_run", "Compute ITV", class = "btn-primary w-100", icon = icon("arrows-left-right")),
        tags$hr(),
        numericInput("acc_perm", "Resampling permutations (accumulation curves)", value = 49, min = 19, step = 10),
        actionButton("acc_run", "Accumulation curves", class = "btn-outline-primary w-100"),
        tags$hr(),
        downloadButton("dl_itv_png", "ITV figure (.png)", class = "w-100 mb-1"),
        downloadButton("dl_acc_png", "Accumulation figure (.png)", class = "w-100 mb-1"),
        downloadButton("dl_itv_csv", "ITV table (.csv)", class = "w-100")
      ),
      step_card("Inter- vs intraspecific variance per trait", icon = "chart-bar",
        plotOutput("itv_plot", height = "480px"),
        verbatimTextOutput("itv_print")
      ),
      step_card("Accumulation of within-species variance", icon = "chart-line",
        plotOutput("acc_plot", height = "480px"),
        verbatimTextOutput("acc_print")
      ),
      step_card("Coefficient of variation per species and trait", icon = "table",
        div(style = "overflow-x:auto;", tableOutput("cv_table"))
      )
    )
  ),

  # ---------------------------------------------------------------- 8. Shape
  nav_panel("8 · Shape", icon = icon("draw-polygon"),
    layout_sidebar(width = sidebar_w,
      sidebar = sidebar(
        tags$h5("Shape space (geometric morphometrics)"),
        help_text("Generalised Procrustes Analysis on the 19 anatomical landmarks (gpa_fish, built on geomorph), then a PCA of the aligned shapes (shape_space)."),
        selectInput("ss_style", "Display", choices = c("Spider + 95% ellipse" = "spider", "Convex hull" = "hull", "Density contour" = "density", "Points only" = "none")),
        checkboxInput("ss_rm_outliers", "Remove Procrustes-distance outliers", FALSE),
        actionButton("ss_run", "Build the shape space", class = "btn-primary w-100", icon = icon("draw-polygon")),
        tags$hr(),
        actionButton("sd_run", "Shape disparity per species (permutation test)", class = "btn-outline-primary w-100"),
        tags$hr(),
        downloadButton("dl_ss_png", "Figure (.png)", class = "w-100 mb-1"),
        downloadButton("dl_ss_pdf", "Figure (.pdf)", class = "w-100")
      ),
      step_card("Shape space", icon = "draw-polygon",
        plotOutput("ss_plot", height = "560px"),
        verbatimTextOutput("ss_print")
      ),
      step_card("Shape disparity", icon = "scale-balanced",
        verbatimTextOutput("sd_print")
      )
    )
  ),

  # ---------------------------------------------------------------- 9. Repeatability
  nav_panel("9 · Repeatability", icon = icon("repeat"),
    layout_sidebar(width = sidebar_w,
      sidebar = sidebar(
        tags$h5("Measurement error"),
        help_text("Needs the same fish digitized several times: the 'bias' sheet of the digitizer workbook (mode = 'repeat'), or the example repeat trial. Percent measurement error and repeatability R (Bailey & Byrnes 1990) per segment, and the placement error of each landmark (digitization_error)."),
        actionButton("rep_example", "Use the example repeat trial", class = "btn-outline-primary w-100 mb-2"),
        fileInput("rep_xlsx", "...or upload a workbook with a 'bias' sheet", accept = c(".xlsx", ".xls"), buttonLabel = "Browse...", placeholder = "landmarks.xlsx"),
        actionButton("rep_run", "Compute repeatability", class = "btn-primary w-100", icon = icon("repeat")),
        tags$hr(),
        downloadButton("dl_rep_png", "Landmark error figure (.png)", class = "w-100 mb-1"),
        downloadButton("dl_rep_csv", "Repeatability table (.csv)", class = "w-100")
      ),
      layout_columns(col_widths = c(5, 7),
        step_card("Repeatability per segment", icon = "table",
          verbatimTextOutput("rep_status"),
          tableOutput("rep_table")
        ),
        step_card("Placement error per landmark", icon = "chart-bar",
          plotOutput("rep_plot", height = "460px"),
          verbatimTextOutput("rep_print")
        )
      )
    )
  ),

  # ---------------------------------------------------------------- 10. Export
  nav_panel("10 · Export", icon = icon("download"),
    layout_sidebar(width = sidebar_w,
      sidebar = sidebar(
        tags$h5("Take everything with you"),
        help_text("A zip with every table (.csv), every figure (.png and .pdf) produced so far, the corrected landmarks, and the R script equivalent to what you clicked."),
        textInput("export_name", "Name for the archive", value = "intraitR_results"),
        downloadButton("dl_all", "Download the results (.zip)", class = "btn-primary w-100")
      ),
      step_card("The R script equivalent to this session", icon = "code",
        help_text("Copy it into RStudio to reproduce (and extend) the analysis in R."),
        verbatimTextOutput("script_view")
      )
    )
  ),

  nav_spacer()
)

# One application, two pages: the pipeline, and -- under ?page=digitizer&job= --
# the package's digitizer for that job. The UI is a function of the request so
# that the digitizer's own page (its sidebar, its theme, its full height) is
# served as it was written, rather than squeezed inside a tab.
ui <- function(req) {
  tok <- job_token_of(req$QUERY_STRING)
  if (!is.null(tok)) {
    e <- job_env(tok)
    if (!is.null(e)) return(e$ui)
    return(fluidPage(tags$h4("This digitizing session has expired."),
                     tags$p("Go back to the pipeline tab, upload the photographs again and open the digitizer.")))
  }
  pipeline_ui
}

# ----------------------------------------------------------------------------- SERVER
server <- function(input, output, session) {

  # ---- the digitizer page: hand the session over to the package's own server
  tok <- job_token_of(isolate(session$clientData$url_search))
  if (!is.null(tok)) {
    e <- job_env(tok)
    if (!is.null(e)) e$server(input, output, session)
    return(invisible(NULL))
  }

  rv <- reactiveValues(
    lm = NULL,          # raw landmarks (intrait_landmarks)
    lm_clean = NULL,    # after corrections
    source = NULL,      # text describing the origin
    photo_dir = NULL,   # folder of photographs (manual digitizer), for backgrounds
    segments = NULL, ratios = NULL, summary = NULL,
    ts = NULL, td = NULL, pf = NULL, itv = NULL, acc = NULL, iv = NULL,
    gpa = NULL, ss = NULL, sd = NULL, outl = NULL,
    rep_lm = NULL, rep_table = NULL, rep_de = NULL, rep_source = NULL,
    code = c("library(intraitR)", ""),
    # digitizing job (the package's own digitizer, opened in a second tab)
    job = NULL, job_base = NULL, digit_refresh = NULL
  )

  add_code <- function(...) rv$code <- c(rv$code, ...)

  active_lm <- reactive(rv$lm_clean %||% rv$lm)

  # ----------------------------------------------------------------- 1. Start
  output$status_box <- renderUI({
    lm <- active_lm()
    if (is.null(lm)) return(tags$p(class = "status-no", icon("circle-xmark"), " No data loaded yet."))
    sp <- species_of(lm)
    tags$div(
      tags$p(class = "status-ok", icon("circle-check"), " ", rv$source),
      tags$ul(
        tags$li(sprintf("%d specimens, %d landmarks", dim(lm$coords)[3], dim(lm$coords)[1])),
        tags$li(if (is.null(sp)) "No species information: group-level analyses will use a single group."
                else sprintf("%d species: %s", length(unique(sp)), paste(sort(unique(sp)), collapse = ", "))),
        tags$li(if (is.null(rv$lm_clean)) "Corrections: not applied yet (tab 3)." else "Corrections applied (tab 3)."),
        tags$li(if (is.null(rv$ratios)) "Traits: not computed yet (tab 4)." else sprintf("Traits computed for %d specimens (tab 4).", nrow(rv$ratios)))
      )
    )
  })

  reset_downstream <- function() {
    rv$lm_clean <- NULL; rv$segments <- NULL; rv$ratios <- NULL; rv$summary <- NULL
    rv$ts <- NULL; rv$td <- NULL; rv$pf <- NULL; rv$itv <- NULL; rv$acc <- NULL; rv$iv <- NULL
    rv$gpa <- NULL; rv$ss <- NULL; rv$sd <- NULL; rv$outl <- NULL
  }

  observeEvent(input$load_example, {
    lm <- safely(load_t26_saudrune_landmarks(), "load_t26_saudrune_landmarks", "Loading the example data...")
    req(lm)
    reset_downstream()
    rv$lm <- lm; rv$photo_dir <- NULL
    rv$source <- "Example data: T-26 La Saudrune survey (intraitR)"
    rv$code <- c("library(intraitR)", "", "## 1. Data", "fish <- load_t26_saudrune_landmarks()")
    showNotification("Example data loaded. Next: tab 3 (Check) or tab 4 (Traits).", type = "message")
  })

  output$sheet_ui <- renderUI({
    req(input$upload_xlsx)
    sheets <- tryCatch(readxl::excel_sheets(input$upload_xlsx$datapath), error = function(e) character(0))
    selectInput("sheet", "Sheet", choices = sheets,
                selected = if ("measurements" %in% sheets) "measurements" else sheets[1])
  })

  observeEvent(input$load_xlsx, {
    req(input$upload_xlsx, input$sheet)
    f <- input$upload_xlsx$datapath
    head <- tryCatch(readxl::read_excel(f, sheet = input$sheet, n_max = 1), error = function(e) NULL)
    if (is.null(head)) { showNotification("Cannot read this sheet.", type = "error"); return() }
    lay <- detect_layout(names(head))
    if (is.null(lay)) {
      showNotification("No landmark columns found: expected 1_X, 1_Y, ... (digitizer) or X_1, Y_1, ...", type = "error"); return()
    }
    idc <- trimws(input$id_col)
    if (!idc %in% names(head)) {
      showNotification(sprintf("Column '%s' not found in the sheet. Available: %s", idc,
                               paste(head(names(head), 12), collapse = ", ")), type = "error"); return()
    }
    lm <- safely(read_landmarks_xlsx(f, sheet = input$sheet, n_landmarks = lay$n,
                                     x_pattern = lay$x, y_pattern = lay$y, id_cols = idc),
                 "read_landmarks_xlsx", "Reading the workbook...")
    req(lm)
    # species from a determinations table (uid, species)
    det_note <- ""
    if (!is.null(input$upload_det)) {
      det <- tryCatch(utils::read.csv(input$upload_det$datapath, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(det) && all(c("uid", "species") %in% names(det))) {
        uid <- sub("_i[0-9]+$", "", lm$metadata$specimen)
        sp  <- det$species[match(uid, det$uid)]
        lm$metadata$species <- sp
        det_note <- sprintf(" + species for %d/%d specimens from the determinations table", sum(!is.na(sp)), length(sp))
      } else showNotification("Determinations table ignored: needs columns 'uid' and 'species'.", type = "warning")
    }
    reset_downstream()
    rv$lm <- lm; rv$photo_dir <- NULL
    rv$source <- paste0("Workbook ", input$upload_xlsx$name, " (sheet ", input$sheet, ")", det_note)
    rv$code <- c("library(intraitR)", "", "## 1. Data",
                 sprintf('fish <- read_landmarks_xlsx("%s", sheet = "%s", n_landmarks = %d,', input$upload_xlsx$name, input$sheet, lay$n),
                 sprintf('                            x_pattern = "%s", y_pattern = "%s", id_cols = "%s")', lay$x, lay$y, idc))
    if (nzchar(det_note)) add_code(sprintf('det <- read.csv("%s")', input$upload_det$name),
                                  'uid <- sub("_i[0-9]+$", "", fish$metadata$specimen)',
                                  'fish$metadata$species <- det$species[match(uid, det$uid)]')
    showNotification(sprintf("Workbook loaded: %d specimens.", dim(lm$coords)[3]), type = "message")
  })

  observeEvent(input$goto_measure, nav_select("nav", "2 · Measure"))

  # ----------------------------------------------------------------- 2. Measure
  # The job: a folder of photographs, a workbook and a journal, handed to the
  # package's digitizer through the `intraitR.digitizer` option (see job_env()).
  observeEvent(input$photos, {
    ph <- input$photos
    token <- paste0("job_", as.integer(Sys.time()), "_", paste(sample(c(letters, 0:9), 6, TRUE), collapse = ""))
    base <- file.path(tempdir(), token)
    dir.create(file.path(base, "photos"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(base, "measurements"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(base, "journal"), recursive = TRUE, showWarnings = FALSE)
    paths <- file.path(base, "photos", ph$name)
    file.copy(ph$datapath, paths, overwrite = TRUE)
    rv$job <- token; rv$job_base <- base
    rv$photo_dir <- file.path(base, "photos")
    showNotification(sprintf("%d photograph(s) ready. Open the digitizer when your settings are set.", nrow(ph)),
                     type = "message")
  })

  # Re-registered whenever a setting changes: the environment is only built when
  # the digitizer page is first opened, so a tab already open keeps its own.
  digit_cfg <- reactive({
    req(rv$job_base)
    photos <- sort(list.files(file.path(rv$job_base, "photos"), full.names = TRUE,
                              pattern = "\\.(jpe?g|png|gif|bmp|tiff?)$", ignore.case = TRUE))
    op <- if (nzchar(trimws(input$operator %||% ""))) trimws(input$operator) else "operator"
    list(photo_dir = file.path(rv$job_base, "photos"), photos = photos,
         xlsx_path = file.path(rv$job_base, "measurements", paste0(op, "_landmarks.xlsx")),
         journal_dir = file.path(rv$job_base, "journal"),
         operator = op, mode = input$digit_mode %||% "new", n_repeats = 3L,
         ruler_mm = max(0.1, as.numeric(input$ruler_mm %||% 10)),
         individuals_per_photo = max(1L, as.integer(input$ind_per_photo %||% 1L)),
         individual_order = "top", individuals_per_row = NA_integer_,
         xlsx_flush_every = 1L, provenance_keys = character(0),
         sheet_measurements = "measurements", sheet_bias = "bias",
         sheet_summary = "bias_summary",
         app_version = paste0("intraitR-pipeline ", APP_VERSION))
  })

  output$digit_open_ui <- renderUI({
    if (is.null(rv$job_base)) return(help_text("Upload photographs first."))
    job_new(rv$job, digit_cfg())      # (re)register with the current settings
    tags$a(href = paste0("?page=digitizer&job=", rv$job), target = "_blank",
           class = "btn btn-primary w-100",
           icon("up-right-from-square"), " Open the digitizer")
  })

  # What the journal holds right now (the workbook is only an export of it).
  digit_journal <- reactive({
    rv$digit_refresh
    req(rv$job_base)
    jd <- file.path(rv$job_base, "journal")
    if (!length(list.files(jd, pattern = "^landmarks_.*\\.tsv$"))) return(NULL)
    tryCatch(consolidate_landmarks(jd, points = 1:N_LM_READ), error = function(e) NULL)
  })

  # Poll the journal while the digitizer is open in the other tab.
  observe({
    req(rv$job_base)
    invalidateLater(4000, session)
    rv$digit_refresh <- Sys.time()
  })

  output$digit_status <- renderPrint({
    if (is.null(rv$job_base)) { cat("No photographs uploaded yet.\n"); return() }
    n_ph <- length(list.files(file.path(rv$job_base, "photos")))
    j <- digit_journal()
    cat(sprintf("Photographs: %d\n", n_ph))
    cat(sprintf("Operator:    %s\n", digit_cfg()$operator))
    if (is.null(j)) { cat("Saved so far: nothing yet (the digitizer writes as you press 'Save & next').\n"); return() }
    sheet <- j$target_sheet; sheet[is.na(sheet) | !nzchar(sheet)] <- "measurements"
    cat(sprintf("Saved so far: %d specimen(s)%s\n", sum(sheet == "measurements"),
                if (any(sheet == "bias")) sprintf(", plus %d repeat digitization(s)", sum(sheet == "bias")) else ""))
  })

  output$digit_progress <- renderTable({
    j <- digit_journal(); req(j)
    sheet <- j$target_sheet; sheet[is.na(sheet) | !nzchar(sheet)] <- "measurements"
    d <- j[sheet == "measurements", , drop = FALSE]
    req(nrow(d))
    data.frame(specimen = d$specimen,
               points_clicked = d$n_clicked,
               never_checked = d$n_seeded,
               unmeasurable = d$n_na,
               scale = ifelse(is.finite(d$mm_per_px), "yes", "MISSING"),
               stringsAsFactors = FALSE)
  }, striped = TRUE, spacing = "xs")

  # Journal -> workbook -> intrait_landmarks, i.e. exactly the route a campaign
  # takes (consolidate_landmarks() then read_landmarks_xlsx()).
  observeEvent(input$digit_load, {
    req(rv$job_base)
    jd <- file.path(rv$job_base, "journal")
    if (!length(list.files(jd, pattern = "^landmarks_.*\\.tsv$"))) {
      showNotification("Nothing saved yet: measure at least one specimen in the digitizer, then come back.", type = "warning"); return()
    }
    res <- safely({
      xl <- file.path(rv$job_base, "measurements", "consolidated.xlsx")
      consolidate_landmarks(jd, points = 1:N_LM_READ, xlsx_path = xl)
      sheets <- readxl::excel_sheets(xl)
      lm <- read_landmarks_xlsx(xl, sheet = "measurements", n_landmarks = N_LM_READ,
                                x_pattern = "{i}_X", y_pattern = "{i}_Y", id_cols = "specimen")
      rep <- if ("bias" %in% sheets)
        tryCatch(read_landmarks_xlsx(xl, sheet = "bias", n_landmarks = N_LM_READ,
                                     x_pattern = "{i}_X", y_pattern = "{i}_Y",
                                     id_cols = c("individual", "operator", "replicate")),
                 error = function(e) NULL) else NULL
      list(lm = lm, rep = rep, xl = xl)
    }, "consolidate_landmarks", "Reading the journal...")
    req(res)
    lm <- res$lm
    # Species are never read from a file name: they come from a determinations
    # table when there is one (uid, species), and are left unknown otherwise.
    det_note <- ""
    if (!is.null(input$digit_det)) {
      det <- tryCatch(utils::read.csv(input$digit_det$datapath, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(det) && all(c("uid", "species") %in% names(det))) {
        uid <- sub("_i[0-9]+$", "", lm$metadata$specimen)
        lm$metadata$species <- det$species[match(uid, det$uid)]
        det_note <- sprintf(", species for %d/%d from the determinations table",
                            sum(!is.na(lm$metadata$species)), nrow(lm$metadata))
      } else showNotification("Determinations table ignored: needs columns 'uid' and 'species'.", type = "warning")
    }
    reset_downstream()
    rv$lm <- lm
    rv$source <- sprintf("Digitized in this session: %d specimen(s)%s", dim(lm$coords)[3], det_note)
    if (!is.null(res$rep)) { rv$rep_lm <- res$rep; rv$rep_source <- "Repeats digitized in this session" }
    rv$code <- c("library(intraitR)", "", "## 1. Data -- digitized with intraitR's own digitizer",
                 '# digitize_landmarks(photo_dir = "photos", xlsx_path = "measurements/<OP>_landmarks.xlsx",',
                 '#                    journal_dir = "journal", operator = "<OP>", mode = "new")',
                 'consolidate_landmarks("journal", points = 1:23, xlsx_path = "measurements/consolidated.xlsx")',
                 'fish <- read_landmarks_xlsx("measurements/consolidated.xlsx", sheet = "measurements", n_landmarks = 23,',
                 '                            x_pattern = "{i}_X", y_pattern = "{i}_Y", id_cols = "specimen")')
    if (nzchar(det_note))
      add_code('det <- read.csv("determinations.csv")',
               'uid <- sub("_i[0-9]+$", "", fish$metadata$specimen)',
               'fish$metadata$species <- det$species[match(uid, det$uid)]')
    showNotification(sprintf("%d specimen(s) loaded. Next: tab 3 (Check).", dim(lm$coords)[3]), type = "message")
  })

  output$digit_download <- downloadHandler(
    filename = function() paste0("intraitR_measurements_", format(Sys.Date(), "%Y%m%d"), ".zip"),
    content = function(file) {
      req(rv$job_base)
      # the workbook is an export; rebuild it from the journal before shipping
      try(consolidate_landmarks(file.path(rv$job_base, "journal"), points = 1:N_LM_READ,
                                xlsx_path = digit_cfg()$xlsx_path), silent = TRUE)
      out <- file.path(tempdir(), paste0("send_", as.integer(Sys.time())))
      dir.create(out, showWarnings = FALSE)
      file.copy(file.path(rv$job_base, "measurements"), out, recursive = TRUE)
      file.copy(file.path(rv$job_base, "journal"), out, recursive = TRUE)
      unlink(file.path(out, "measurements", "consolidated.xlsx"))
      unlink(list.files(file.path(out, "measurements"), pattern = "\\.prev\\.xlsx$", full.names = TRUE))
      writeLines(c("Measurements made with the intraitR pipeline app.",
                   paste("Operator:", digit_cfg()$operator),
                   paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M")),
                   paste("intraitR", as.character(utils::packageVersion("intraitR"))),
                   "",
                   "measurements/ : the workbook (an export).",
                   "journal/      : the append-only journal -- the source of truth.",
                   "Send BOTH; the journal is what gets merged with consolidate_landmarks()."),
                 file.path(out, "README.txt"))
      zip_dir(out, file)
    }
  )

  zip_dir <- function(dir, file) {
    old <- setwd(dir); on.exit(setwd(old))
    files <- list.files(".", recursive = TRUE)
    if (has_pkg("zip")) zip::zip(file, files) else utils::zip(file, files, flags = "-r9Xq")
  }

  # ----------------------------------------------------------------- 3. Check
  output$qc_specimen_ui <- renderUI({
    lm <- rv$lm; req(lm)
    selectInput("qc_spec", "Specimen", choices = lm$metadata$specimen %||% dimnames(lm$coords)[[3]])
  })

  output$qc_plot <- renderPlot({
    lm <- rv$lm; req(lm, input$qc_spec)
    bg <- if (isTRUE(input$qc_photo)) photo_for_code(rv$photo_dir, input$qc_spec) else NULL
    obj <- if (is.null(bg) && !is.null(rv$lm_clean)) rv$lm_clean else lm
    tryCatch(plot_fishmorph_points(obj, specimen = input$qc_spec, background_image = bg,
                                   main = if (is.null(bg) && !is.null(rv$lm_clean)) "corrected coordinates" else "raw coordinates"),
             error = function(e) { plot.new(); text(0.5, 0.5, conditionMessage(e)) })
  })

  observeEvent(input$qc_apply, {
    lm <- rv$lm; req(lm)
    code <- c("", "## 3. Quality control")
    x <- lm
    withProgress(message = "Applying corrections...", value = 0.2, {
      if (isTRUE(input$qc_impute)) { x <- tryCatch(impute_landmarks(x), error = function(e) { notify_error(e, "impute_landmarks"); x }); code <- c(code, "fish <- impute_landmarks(fish)") }
      incProgress(0.3)
      if (isTRUE(input$qc_orient)) { x <- tryCatch(standardize_orientation(x), error = function(e) { notify_error(e, "standardize_orientation"); x }); code <- c(code, "fish <- standardize_orientation(fish)") }
      incProgress(0.3)
      if (isTRUE(input$qc_geom))   { x <- tryCatch(correct_geometry(x), error = function(e) { notify_error(e, "correct_geometry"); x }); code <- c(code, "fish <- correct_geometry(fish)") }
    })
    rv$lm_clean <- x
    rv$segments <- NULL; rv$ratios <- NULL
    add_code(code)
    showNotification("Corrections applied. Next: tab 4 (Traits).", type = "message")
  })

  observeEvent(input$qc_outliers, {
    lm <- active_lm(); req(lm)
    res <- safely({
      sh <- fishmorph_shape_landmarks(lm, drop_incomplete = TRUE)
      g  <- gpa_fish(sh, flag_outliers = FALSE)
      detect_outliers(g, threshold = input$qc_thr, plot = FALSE)
    }, "detect_outliers", "Procrustes superimposition and outlier screening...")
    req(res)
    rv$outl <- res
    add_code("", sprintf("outl <- detect_outliers(gpa_fish(fishmorph_shape_landmarks(fish)), threshold = %s, plot = FALSE)", input$qc_thr))
  })

  output$qc_summary <- renderPrint({
    lm <- rv$lm_clean
    if (is.null(lm)) { cat("No correction applied yet.\n"); return() }
    co <- lm$coords
    imp <- attr(co, "imputed"); cor <- attr(co, "corrected"); ori <- attr(co, "orientation_log")
    cat(sprintf("Imputed points:   %s\n", if (is.null(imp)) "none" else sum(imp, na.rm = TRUE)))
    cat(sprintf("Corrected points: %s\n", if (is.null(cor)) "none" else sum(cor, na.rm = TRUE)))
    if (!is.null(ori)) cat(sprintf("Mirrored specimens: %d horizontally, %d vertically\n", sum(ori$flipped_x), sum(ori$flipped_y)))
    if (!is.null(rv$outl)) cat(sprintf("\nShape outliers (threshold %s): %d flagged\n", input$qc_thr, length(rv$outl$outliers)))
  })

  output$qc_outlier_table <- renderTable({
    req(rv$outl)
    r <- rv$outl$rank; r <- r[r$outlier %in% TRUE, , drop = FALSE]
    if (!nrow(r)) return(data.frame(note = "No outlier flagged."))
    r$procrustes_distance <- round(r$procrustes_distance, 4)
    head(r[, c("specimen", "procrustes_distance")], 30)
  }, striped = TRUE, spacing = "xs")

  # ----------------------------------------------------------------- 4. Traits
  observeEvent(input$traits_run, {
    lm <- active_lm(); req(lm)
    sp <- species_of(lm)
    seg <- safely(fishmorph_segments(lm, scale_cm = input$scale_cm, scale_action = input$scale_action),
                  "fishmorph_segments", "Computing the 11 segments...")
    req(seg)
    rat <- safely(fishmorph_ratios(seg, no_caudal_fin = input$no_caudal, ventral_mouth = input$ventral_mouth,
                                   no_pectoral_fin = input$no_pectoral, na_action = input$na_action_tr,
                                   groups = if (input$na_action_tr == "impute_group_mean") sp else NULL),
                  "fishmorph_ratios", "Computing the 9 ratios...")
    req(rat)
    rv$segments <- seg; rv$ratios <- rat
    spr <- if (!is.null(rat$species)) rat$species else if (!is.null(sp)) sp[match(rownames(rat), lm$metadata$specimen)] else NULL
    rv$summary <- if (!is.null(spr) && !all(is.na(spr))) tryCatch(summary_traits(rat[, RATIO_COLS], groups = spr), error = function(e) NULL) else NULL
    rv$ts <- NULL; rv$td <- NULL; rv$pf <- NULL; rv$itv <- NULL; rv$acc <- NULL; rv$iv <- NULL
    add_code("", "## 4. FISHMORPH traits",
             sprintf('segments <- fishmorph_segments(fish, scale_cm = %s, scale_action = "%s")', input$scale_cm, input$scale_action),
             sprintf('ratios   <- fishmorph_ratios(segments, no_caudal_fin = %s, ventral_mouth = %s, no_pectoral_fin = %s, na_action = "%s")',
                     input$no_caudal, input$ventral_mouth, input$no_pectoral, input$na_action_tr),
             'ratio_cols <- c("BEl", "VEp", "REs", "OGp", "RMl", "BLs", "PFv", "PFs", "CPt")',
             'ok      <- complete.cases(ratios[ratio_cols])',
             'traits  <- ratios[ok, ratio_cols]',
             'species <- ratios$species[ok]',
             'summary_traits(traits, groups = species)')
    showNotification(sprintf("Traits computed for %d specimens. Next: tabs 5-8.", nrow(rat)), type = "message")
  })

  # traits ready for the ordinations: complete cases of the 9 ratios + species
  traits_ready <- reactive({
    req(rv$ratios)
    r <- rv$ratios
    ok <- stats::complete.cases(r[, RATIO_COLS])
    sp <- if (!is.null(r$species)) r$species else NULL
    if (is.null(sp) || all(is.na(sp))) sp <- rep("all specimens", nrow(r))
    sp[is.na(sp)] <- "unknown"
    list(x = r[ok, RATIO_COLS], groups = sp[ok], n_dropped = sum(!ok))
  })

  output$traits_msg <- renderPrint({
    if (is.null(rv$ratios)) { cat("Click 'Compute the traits'.\n"); return() }
    tr <- traits_ready()
    cat(sprintf("%d specimens; %d with at least one missing ratio (excluded from the ordinations); %d groups.\n",
                nrow(rv$ratios), tr$n_dropped, length(unique(tr$groups))))
  })
  output$ratios_table <- renderTable({
    req(rv$ratios); r <- rv$ratios
    keep <- intersect(c("specimen", "species", RATIO_COLS), names(r))
    head(ratio_table(r[, keep, drop = FALSE]), 200)
  }, striped = TRUE, spacing = "xs")
  output$summary_table <- renderTable({ req(rv$summary); ratio_table(rv$summary) }, striped = TRUE, spacing = "xs")

  output$dl_segments <- downloadHandler("fishmorph_segments.csv", function(f) utils::write.csv(rv$segments, f, row.names = FALSE))
  output$dl_ratios   <- downloadHandler("fishmorph_ratios.csv",   function(f) utils::write.csv(rv$ratios, f, row.names = FALSE))
  output$dl_summary  <- downloadHandler("traits_summary_by_species.csv", function(f) utils::write.csv(rv$summary, f, row.names = FALSE))

  # ----------------------------------------------------------------- 5. Trait space
  observeEvent(input$ts_run, {
    tr <- traits_ready()
    ax <- as.integer(strsplit(input$ts_axes, ",")[[1]])
    ts <- safely(trait_space(tr$x, groups = tr$groups, axes = ax, remove_outliers = input$ts_outliers),
                 "trait_space", "Building the trait space...")
    req(ts); rv$ts <- ts
    add_code("", "## 5. Trait space",
             sprintf("ts <- trait_space(traits, groups = species, axes = c(%d, %d), remove_outliers = %s)", ax[1], ax[2], input$ts_outliers),
             sprintf('plot(ts, style = "%s", legend_italic = %s)', input$ts_style, input$ts_italic))
  })
  draw_ts <- function() { req(rv$ts); plot(rv$ts, style = input$ts_style, legend_italic = input$ts_italic) }
  output$ts_plot  <- renderPlot(draw_ts())
  output$ts_print <- renderPrint({ req(rv$ts); print(rv$ts) })

  observeEvent(input$td_run, {
    tr <- traits_ready()
    td <- safely(trait_disparity(tr$x, groups = tr$groups, iter = input$td_iter), "trait_disparity", "Permutation test...")
    req(td); rv$td <- td
    add_code(sprintf("td <- trait_disparity(traits, groups = species, iter = %d)", input$td_iter))
  })
  output$td_print <- renderPrint({ if (is.null(rv$td)) cat("Not computed yet.\n") else print(rv$td) })

  fig_download <- function(draw, name) {
    list(
      png = downloadHandler(paste0(name, ".png"), function(f) { png(f, 2000, 1500, res = 220); draw(); dev.off() }),
      pdf = downloadHandler(paste0(name, ".pdf"), function(f) { pdf(f, 9, 7); draw(); dev.off() })
    )
  }
  d <- fig_download(draw_ts, "trait_space"); output$dl_ts_png <- d$png; output$dl_ts_pdf <- d$pdf

  # ----------------------------------------------------------------- 6. FISHMORPH
  observeEvent(input$pf_run, {
    tr <- traits_ready()
    pf <- safely(project_fishmorph(tr$x, groups = tr$groups), "project_fishmorph", "Projecting into the FISHMORPH space (~9 000 species)...")
    req(pf); rv$pf <- pf
    add_code("", "## 6. FISHMORPH morphospace",
             "proj <- project_fishmorph(traits, groups = species)",
             sprintf('plot(proj, style = "%s", itv_reference = %s, arrows = %s)', input$pf_style, input$pf_itvref, input$pf_arrows))
  })
  draw_pf <- function() { req(rv$pf); plot(rv$pf, style = input$pf_style, itv_reference = input$pf_itvref, arrows = input$pf_arrows) }
  output$pf_plot  <- renderPlot(draw_pf())
  output$pf_print <- renderPrint({ req(rv$pf); print(rv$pf) })
  d <- fig_download(draw_pf, "fishmorph_projection"); output$dl_pf_png <- d$png; output$dl_pf_pdf <- d$pdf

  # ----------------------------------------------------------------- 7. ITV
  observeEvent(input$itv_run, {
    tr <- traits_ready()
    if (length(unique(tr$groups)) < 2) showNotification("ITV needs at least two species.", type = "warning")
    itv <- safely(itv_index(tr$x, groups = tr$groups), "itv_index", "Partitioning the variance...")
    req(itv); rv$itv <- itv
    rv$iv <- tryCatch(intraspecific_variability(traits = tr$x, groups = tr$groups), error = function(e) NULL)
    add_code("", "## 7. Intraspecific variability",
             "itv <- itv_index(traits, groups = species); plot(itv)",
             "iv  <- intraspecific_variability(traits = traits, groups = species)   # CV per species x trait")
  })
  draw_itv <- function() { req(rv$itv); plot(rv$itv) }
  output$itv_plot  <- renderPlot(draw_itv())
  output$itv_print <- renderPrint({ req(rv$itv); print(rv$itv) })
  output$cv_table  <- renderTable({ req(rv$iv); ratio_table(rv$iv$trait_cv) }, striped = TRUE, spacing = "xs")

  observeEvent(input$acc_run, {
    tr <- traits_ready()
    acc <- safely(itv_accumulation(tr$x, groups = tr$groups, n_perm = input$acc_perm), "itv_accumulation", "Resampling (this can take a minute)...")
    req(acc); rv$acc <- acc
    add_code(sprintf("acc <- itv_accumulation(traits, groups = species, n_perm = %d); plot(acc)", input$acc_perm))
  })
  draw_acc <- function() { req(rv$acc); plot(rv$acc) }
  output$acc_plot  <- renderPlot(draw_acc())
  output$acc_print <- renderPrint({ if (is.null(rv$acc)) cat("Not computed yet.\n") else print(rv$acc) })
  output$dl_itv_png <- downloadHandler("itv_index.png", function(f) { png(f, 2000, 1500, res = 220); draw_itv(); dev.off() })
  output$dl_acc_png <- downloadHandler("itv_accumulation.png", function(f) { png(f, 2000, 1500, res = 220); draw_acc(); dev.off() })
  output$dl_itv_csv <- downloadHandler("itv_index.csv", function(f) utils::write.csv(rv$itv$per_trait, f, row.names = FALSE))

  # ----------------------------------------------------------------- 8. Shape
  observeEvent(input$ss_run, {
    lm <- active_lm(); req(lm)
    res <- safely({
      sh <- fishmorph_shape_landmarks(lm, drop_incomplete = TRUE)
      g  <- gpa_fish(sh, remove_outliers = input$ss_rm_outliers)
      sp <- if (!is.null(g$metadata$species)) g$metadata$species else species_of(sh)
      if (is.null(sp) || all(is.na(sp))) sp <- NULL
      list(gpa = g, ss = shape_space(g, groups = sp))
    }, "shape_space", "Procrustes superimposition and PCA...")
    req(res); rv$gpa <- res$gpa; rv$ss <- res$ss; rv$sd <- NULL
    add_code("", "## 8. Shape space",
             "shape <- fishmorph_shape_landmarks(fish)",
             sprintf("gpa   <- gpa_fish(shape, remove_outliers = %s)", input$ss_rm_outliers),
             "ss    <- shape_space(gpa, groups = gpa$metadata$species)",
             sprintf('plot(ss, style = "%s")', input$ss_style))
  })
  draw_ss <- function() { req(rv$ss); plot(rv$ss, style = input$ss_style) }
  output$ss_plot  <- renderPlot(draw_ss())
  output$ss_print <- renderPrint({ req(rv$ss); print(rv$ss); cat("\n"); print(rv$gpa) })
  d <- fig_download(draw_ss, "shape_space"); output$dl_ss_png <- d$png; output$dl_ss_pdf <- d$pdf

  observeEvent(input$sd_run, {
    req(rv$gpa)
    sp <- rv$gpa$metadata$species
    if (is.null(sp) || length(unique(sp)) < 2) { showNotification("Shape disparity needs at least two species.", type = "warning"); return() }
    sd <- safely(intraspecific_variability(gpa = rv$gpa, groups = sp, iter = 199), "intraspecific_variability", "Permutation test on shape disparity...")
    req(sd); rv$sd <- sd
    add_code("sdisp <- intraspecific_variability(gpa = gpa, groups = gpa$metadata$species, iter = 199)")
  })
  output$sd_print <- renderPrint({ if (is.null(rv$sd)) cat("Not computed yet.\n") else print(rv$sd) })

  # ----------------------------------------------------------------- 9. Repeatability
  observeEvent(input$rep_example, {
    lm <- safely(load_t26_saudrune_landmarks("repeatability"), "load_t26_saudrune_landmarks", "Loading the repeat trial...")
    req(lm); rv$rep_lm <- lm; rv$rep_source <- "Example repeat trial (T-26)"; rv$rep_table <- NULL; rv$rep_de <- NULL
  })
  observeEvent(input$rep_xlsx, {
    f <- input$rep_xlsx$datapath
    sheets <- tryCatch(readxl::excel_sheets(f), error = function(e) character(0))
    if (!"bias" %in% sheets) { showNotification("No 'bias' sheet in this workbook.", type = "error"); return() }
    lm <- safely(read_landmarks_xlsx(f, sheet = "bias", n_landmarks = N_LM_READ, x_pattern = "{i}_X", y_pattern = "{i}_Y",
                                     id_cols = c("individual", "operator", "replicate")),
                 "read_landmarks_xlsx", "Reading the bias sheet...")
    req(lm); rv$rep_lm <- lm; rv$rep_source <- paste0("Sheet 'bias' of ", input$rep_xlsx$name); rv$rep_table <- NULL; rv$rep_de <- NULL
  })

  observeEvent(input$rep_run, {
    lm <- rv$rep_lm; req(lm)
    res <- safely({
      seg <- fishmorph_segments(lm, scale_cm = input$scale_cm)   # scale_action = "na": a replicate without scale bar is dropped, never mixed in pixels
      ind <- lm$metadata$individual[match(rownames(seg), lm$metadata$specimen)]
      tab <- do.call(rbind, lapply(SEG_COLS, function(s) {
        ok <- is.finite(seg[[s]])
        me <- tryCatch(measurement_error(data.frame(individual = ind[ok], value = seg[[s]][ok]), individual = "individual", method = "anova"),
                       error = function(e) NULL)
        data.frame(segment = s, individuals = length(unique(ind[ok])),
                   `%ME` = if (is.null(me)) NA else round(me$percent_measurement_error, 2),
                   R = if (is.null(me)) NA else round(me$repeatability, 3), check.names = FALSE)
      }))
      sp <- if (!is.null(lm$metadata$species)) lm$metadata$species else NULL
      de <- digitization_error(lm, individual = lm$metadata$individual, species = sp, exclude_landmarks = c(20, 21, 25))
      list(tab = tab, de = de)
    }, "measurement_error", "Computing repeatability...")
    req(res); rv$rep_table <- res$tab; rv$rep_de <- res$de
    add_code("", "## 9. Repeatability (repeat trial)",
             if (grepl("Example", rv$rep_source %||% "")) 'rep <- load_t26_saudrune_landmarks("repeatability")'
             else 'rep <- read_landmarks_xlsx("<workbook>.xlsx", sheet = "bias", n_landmarks = 23, x_pattern = "{i}_X", y_pattern = "{i}_Y", id_cols = c("individual", "operator", "replicate"))',
             'seg_rep <- fishmorph_segments(rep)',
             'ind     <- rep$metadata$individual[match(rownames(seg_rep), rep$metadata$specimen)]',
             'measurement_error(data.frame(individual = ind, value = seg_rep$Bl), individual = "individual", method = "anova")',
             'derr <- digitization_error(rep, individual = rep$metadata$individual, exclude_landmarks = c(20, 21, 25)); plot(derr)')
  })
  output$rep_status <- renderPrint({ cat(if (is.null(rv$rep_lm)) "No repeat data loaded.\n" else paste0(rv$rep_source, ": ", dim(rv$rep_lm$coords)[3], " digitizations of ", length(unique(rv$rep_lm$metadata$individual)), " individuals\n")) })
  output$rep_table <- renderTable({ req(rv$rep_table); rv$rep_table }, striped = TRUE, spacing = "xs")
  draw_rep <- function() { req(rv$rep_de); plot(rv$rep_de) }
  output$rep_plot  <- renderPlot(draw_rep())
  output$rep_print <- renderPrint({ req(rv$rep_de); print(rv$rep_de) })
  output$dl_rep_png <- downloadHandler("digitization_error.png", function(f) { png(f, 2000, 1500, res = 220); draw_rep(); dev.off() })
  output$dl_rep_csv <- downloadHandler("repeatability_by_segment.csv", function(f) utils::write.csv(rv$rep_table, f, row.names = FALSE))

  # ----------------------------------------------------------------- 10. Export
  output$script_view <- renderText(paste(rv$code, collapse = "\n"))

  output$dl_all <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_-]", "_", input$export_name), "_", format(Sys.Date(), "%Y%m%d"), ".zip"),
    content = function(file) {
      out <- file.path(tempdir(), paste0("export_", as.integer(Sys.time()))); dir.create(out)
      withProgress(message = "Writing tables, figures and script...", value = 0.2, {
        wcsv <- function(x, n) if (!is.null(x)) utils::write.csv(x, file.path(out, n), row.names = FALSE)
        wcsv(rv$segments, "fishmorph_segments.csv"); wcsv(rv$ratios, "fishmorph_ratios.csv")
        wcsv(rv$summary, "traits_summary_by_species.csv")
        if (!is.null(rv$itv)) { wcsv(rv$itv$per_trait, "itv_per_trait.csv"); wcsv(rv$itv$multivariate, "itv_multivariate.csv") }
        if (!is.null(rv$iv)) wcsv(rv$iv$trait_cv, "cv_by_species_trait.csv")
        if (!is.null(rv$acc)) wcsv(rv$acc$summary, "itv_accumulation_summary.csv")
        if (!is.null(rv$td)) wcsv(data.frame(species = names(rv$td$disparity), disparity = rv$td$disparity), "trait_disparity.csv")
        if (!is.null(rv$ts)) wcsv(cbind(specimen = rownames(rv$ts$scores), rv$ts$scores), "trait_space_scores.csv")
        if (!is.null(rv$ss)) wcsv(cbind(specimen = rownames(rv$ss$scores), rv$ss$scores), "shape_space_scores.csv")
        if (!is.null(rv$pf)) wcsv(cbind(specimen = rownames(rv$pf$scores), rv$pf$scores), "fishmorph_projection_scores.csv")
        wcsv(rv$rep_table, "repeatability_by_segment.csv")
        if (!is.null(rv$outl)) wcsv(rv$outl$rank, "shape_outliers.csv")
        lm <- active_lm()
        if (!is.null(lm)) {
          co <- lm$coords; p <- dim(co)[1]
          wide <- as.data.frame(matrix(aperm(co, c(3, 2, 1)), nrow = dim(co)[3]))
          names(wide) <- as.vector(t(outer(seq_len(p), c("X", "Y"), function(i, a) paste0(i, "_", a))))
          wcsv(cbind(lm$metadata, wide), "landmarks_used.csv")
        }
        incProgress(0.4)
        figs <- list(trait_space = draw_ts, fishmorph_projection = draw_pf, itv_index = draw_itv,
                     itv_accumulation = draw_acc, shape_space = draw_ss, digitization_error = draw_rep)
        for (n in names(figs)) {
          ok <- tryCatch({ png(file.path(out, paste0(n, ".png")), 2000, 1500, res = 220); figs[[n]](); dev.off(); TRUE },
                         error = function(e) { try(dev.off(), silent = TRUE); FALSE })
          if (!ok) unlink(file.path(out, paste0(n, ".png")))
          else { pdf(file.path(out, paste0(n, ".pdf")), 9, 7); try(figs[[n]](), silent = TRUE); dev.off() }
        }
        incProgress(0.3)
        writeLines(c(paste0("## intraitR pipeline app v", APP_VERSION, " -- ", format(Sys.time(), "%Y-%m-%d %H:%M")),
                     paste0("## Data: ", rv$source %||% "none"), "", rv$code), file.path(out, "pipeline_script.R"))
        writeLines(c("intraitR pipeline -- export", paste("Date:", Sys.time()), paste("Data:", rv$source %||% "none"),
                     paste("intraitR", as.character(utils::packageVersion("intraitR"))),
                     paste("Rfishmorph", as.character(utils::packageVersion("Rfishmorph"))),
                     paste("geomorph", as.character(utils::packageVersion("geomorph"))), "",
                     "Files:", list.files(out)), file.path(out, "README.txt"))
      })
      zip_dir(out, file)
    }
  )
}

shinyApp(ui, server)
