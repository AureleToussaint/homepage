# =============================================================================
# DAY 2  ·  GROUP 1  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q1 exercise for GROUP 1.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 3.
#
# RULES
#   • Use the ANNOTATED script (`day2_functional_diversity_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group1_day2_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Compute CWM for the 5 plant traits, link community-level dominant trait values to environmental gradients, and partition the unique effect of disturbance using partial regression.
#
# DATA / OBJECTS YOU NEED
#   • outputs/day1_sp_traits_clean.rds
#   • data/02_environment.csv
#   • data/03_community_composition.csv
#   • trait_mat_sc + comm_mat (built by day2_functional_diversity_annotated.R)
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day2_functional_diversity_annotated.R` (sections relevant to Q1) for context.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. SETUP — load packages and rebuild required objects
# -----------------------------------------------------------------------------
library(tidyverse)
library(FD)
library(funspace)
library(vegan)
library(ggrepel)
library(ggforce)
library(patchwork)
library(corrplot)

theme_set(theme_bw(base_size = 12))
if (!dir.exists("outputs")) dir.create("outputs")

# Pre-requisite: outputs/day1_sp_traits_clean.rds must exist (Day 1 output).
sp_traits <- readRDS("outputs/day1_sp_traits_clean.rds")
env       <- read_csv("data/02_environment.csv",           show_col_types = FALSE)
comm_raw  <- read_csv("data/03_community_composition.csv", show_col_types = FALSE)

# Build matrices (same as in the annotated Day 2 script)
trait_mat <- sp_traits |>
  select(species, SLA, LDMC, LeafN, Height, SeedMass) |>
  mutate(across(c(SLA, LDMC, Height, SeedMass), log)) |>
  column_to_rownames("species") |> as.matrix()
trait_mat_sc <- scale(trait_mat)
comm_mat <- comm_raw |> column_to_rownames("plot_id") |> as.matrix()

# Common FD output (used by Q1-Q4)
fd_full <- dbFD(trait_mat_sc, comm_mat, calc.FRic = TRUE, calc.FDiv = TRUE,
                stand.x = FALSE, m = 3, print.pco = FALSE, messages = FALSE)

# =============================================================================
# Q1 — Community-Weighted Means (CWM) and the mass-ratio hypothesis
# =============================================================================
# -----------------------------------------------------------------------------
# Q1a — CWM(Height) along the disturbance gradient
# -----------------------------------------------------------------------------
# TASK
#   Use `FD::functcomp(trait_mat_sc, comm_mat)` to compute all CWMs. Plot
#   CWM(Height) versus the disturbance index. Colour by elevation. Which
#   PFT(s) drive the pattern at the disturbed end?
#
# EXPECTED OUTPUT
#   Scatter + smooth saved as outputs/group1_day2_Q1a_CWM_Height.png.
#   Disturbance often leads to taller dominants in the recovery range
#   (ruderals/competitors take over).
#
# HINTS
#   • `functcomp()` returns a sites × traits matrix.
#   • Join with env using `bind_cols()` after matching `plot_id`.
#   • Use `geom_smooth(method = 'loess')` to capture non-linearity.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q1b — Correlation matrix CWMs × environment
# -----------------------------------------------------------------------------
# TASK
#   Build a 5 × 5 correlation matrix between the 5 CWMs and 5 environmental
#   variables. Visualise with `corrplot()`. Identify the strongest CWM ×
#   environment association.
#
# EXPECTED OUTPUT
#   PNG outputs/group1_day2_Q1b_corr.png + a sentence stating the strongest
#   pair (e.g. CWM(SLA) × elevation, expected r ≈ -0.6).
#
# HINTS
#   • `cor()` between two matrices/data frames returns the cross matrix.
#   • `corrplot::corrplot(method = 'color', type = 'full', addCoef.col =
#     'black')`.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q1c — Partial regression: CWM(SLA) ~ elevation + disturbance
# -----------------------------------------------------------------------------
# TASK
#   Fit `lm(CWM_SLA ~ elevation + disturbance)`. Test whether disturbance
#   retains a significant effect after controlling for elevation. Discuss in
#   terms of the response–effect framework (which trait responds, which
#   environmental driver dominates).
#
# EXPECTED OUTPUT
#   Model summary with two coefficient p-values; partial regression plots;
#   interpretation paragraph.
#
# HINTS
#   • Use `car::avPlots(mod)` for added-variable plots.
#   • If both predictors are significant, report standardised betas to compare
#     effect sizes.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group1_day2_*.png
# =============================================================================
