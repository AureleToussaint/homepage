# Deploy the intraitR pipeline app to shinyapps.io.
# Run once: rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
# (copied from https://www.shinyapps.io/admin/#/tokens)

# magick is what resizes the photographs on arrival (in C, EXIF rotation
# included). Without it the app falls back to a base-R subsample, so it is
# not required -- but it must be installed HERE for rsconnect to ship it.
pkgs <- c("shiny", "bslib", "geomorph", "readxl", "writexl", "jpeg", "png", "zip",
          "magick", "rsconnect", "remotes")
miss <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(miss)) install.packages(miss)

# intraitR and Rfishmorph must be installed from GitHub so that rsconnect can
# reinstall them on the server (it records the remote origin).
if (!requireNamespace("Rfishmorph", quietly = TRUE)) remotes::install_github("FunTraits/Rfishmorph")
if (!requireNamespace("intraitR",   quietly = TRUE)) remotes::install_github("FunTraits/intraitR")

rsconnect::deployApp(
  appDir   = normalizePath(dirname(if (interactive()) rstudioapi::getActiveDocumentContext()$path else "deploy.R")),
  appName  = "intraitR-pipeline",
  appFiles = c("app.R", "README.md"),
  forceUpdate = TRUE
)
