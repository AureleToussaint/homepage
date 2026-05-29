# =============================================================================
# DAY 2  ·  GROUP 2  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q2 exercise for GROUP 2.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 3.
#
# RULES
#   • Use the ANNOTATED script (`day2_functional_diversity_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group2_day2_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Demonstrate how FRic is confounded with species richness, why FDis is more robust, and how trait number alters FD rankings.
#
# DATA / OBJECTS YOU NEED
#   • fd_full (output of dbFD on the full 5-trait set)
#   • trait_mat_sc, comm_mat
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day2_functional_diversity_annotated.R` (sections relevant to Q2) for context.
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
# Q2 — FD indices: sensitivity to species richness and trait number
# =============================================================================
# -----------------------------------------------------------------------------
# Q2a — Plots with extreme FDis values
# -----------------------------------------------------------------------------
# TASK
#   From `fd_full$FDis`, identify the 3 plots with the highest and the 3 with
#   the lowest FDis. Report their elevation, disturbance and species richness.
#   Do extreme-FDis plots cluster at certain environmental conditions?
#
# EXPECTED OUTPUT
#   Two short tables (top-3 and bottom-3) with their environmental context.
#   Optionally, a small map of plot positions in the elevation × disturbance
#   plane.
#
# HINTS
#   • `tibble(plot = names(fd_full$FDis), FDis = fd_full$FDis) |>
#     left_join(env)`.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q2b — SR confounding: correlation with FRic vs FDis
# -----------------------------------------------------------------------------
# TASK
#   Compute Pearson correlations between species richness (SR) and (i) FRic,
#   (ii) FDis. Which index is more confounded? Propose a method to disentangle
#   FD from SR (residuals or null model).
#
# EXPECTED OUTPUT
#   Two correlation values + a 3-line proposal for SR-correction.
#
# HINTS
#   • Expected: |r(SR, FRic)| > 0.7; |r(SR, FDis)| < 0.3.
#   • Correction options: residuals of FRic ~ SR; null-model SES values (e.g.
#     {picante::ses.mpd}); randomly subsample to fixed S.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q2c — FRic sensitivity to trait number
# -----------------------------------------------------------------------------
# TASK
#   Re-run `dbFD()` keeping only SLA + LDMC (2 traits) and compare to the full
#   5-trait set. Compute Spearman rank correlation between the two FRic
#   vectors. Quantify how much site rankings change.
#
# EXPECTED OUTPUT
#   Spearman r value; scatterplot of FRic_5traits vs FRic_2traits with y=x
#   reference line. Discussion: is the choice of traits ecologically
#   meaningful for the question?
#
# HINTS
#   • `dbFD(trait_mat_sc[, c('SLA','LDMC')], comm_mat, ...)`.
#   • Mouillot et al. (2021) recommend 4–6 trait dimensions and stability
#     checks.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group2_day2_*.png
# =============================================================================
