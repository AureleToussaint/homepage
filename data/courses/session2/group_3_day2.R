# =============================================================================
# DAY 2  ·  GROUP 3  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q3 exercise for GROUP 3.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 3.
#
# RULES
#   • Use the ANNOTATED script (`day2_functional_diversity_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group3_day2_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Move beyond the convex hull. Build a kernel-density TPD with {funspace}, map a third variable (IUCN threat) onto the trait space, and test the observed FRic against a null expectation.
#
# DATA / OBJECTS YOU NEED
#   • trait_mat_sc, sp_traits (with PFT_label and IUCN_status)
#   • funspace package (>= 0.2.0)
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day2_functional_diversity_annotated.R` (sections relevant to Q3) for context.
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
# Q3 — Probabilistic trait spaces with {funspace}
# =============================================================================
# -----------------------------------------------------------------------------
# Q3a — Build the funspace object
# -----------------------------------------------------------------------------
# TASK
#   Run `princomp(trait_mat_sc)` then `funspace(pca, PCs = c(1,2), group.vec =
#   sp_traits$PFT_label, threshold = 0.95)`. Print the summary and plot per
#   group. Which PFT occupies the largest TPD?
#
# EXPECTED OUTPUT
#   summary(fs) printout + outputs/group3_day2_Q3a_funspace.png. Ruderals or
#   competitors typically occupy the largest TPD volume.
#
# HINTS
#   • Make sure `funspace` ≥ 0.2.0 — check with `packageVersion('funspace')`.
#   • `plot(fs, type = 'groups', quant.plot = TRUE, pnt = TRUE)`.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q3b — Map IUCN threat with funspaceGAM
# -----------------------------------------------------------------------------
# TASK
#   Recode IUCN status into Mooers et al. (2008) weights (LC=0, NT=0.1,
#   VU=0.4, EN=0.667, CR=0.999). Run `funspaceGAM(fs, var = iucn_w, family =
#   gaussian())`. Where in the trait space are threatened species
#   concentrated?
#
# EXPECTED OUTPUT
#   GAM map saved as outputs/group3_day2_Q3b_threat.png. Threat should peak at
#   one extreme of LES (typically the conservative corner).
#
# HINTS
#   • Build the weight vector with `recode()` against IUCN_status.
#   • `plot(fs_gam, quant.plot = TRUE)` — red = high predicted threat.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q3c — Null model with funspaceNull
# -----------------------------------------------------------------------------
# TASK
#   Run `funspaceNull(fs, iterations = 999, method = 'mvnorm')`. Is observed
#   FRic significantly different from random? What does a significant test
#   mean ecologically?
#
# EXPECTED OUTPUT
#   Summary with effect size and p-value. Discuss in 4-5 lines: p < 0.05
#   indicates non-random functional structure (e.g. filtering or limiting
#   similarity).
#
# HINTS
#   • If FRic_obs > random → trait space larger than expected → limiting
#     similarity (niche differentiation).
#   • If FRic_obs < random → trait space smaller → environmental filtering.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group3_day2_*.png
# =============================================================================
