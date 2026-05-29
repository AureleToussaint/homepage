# =============================================================================
# DAY 5  ·  GROUP 3  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q3 exercise for GROUP 3.
# Submit your annotated R script with figures saved to outputs/ before Friday evening · final submission.
#
# RULES
#   • Use the ANNOTATED script (`day5_gm_itv_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group3_day5_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Test whether habitat (fast vs slow water) significantly predicts shape, compute pairwise Procrustes distances and visualise residual shape variation after removing the species effect.
#
# DATA / OBJECTS YOU NEED
#   • gpa, pca, meta from Q1
#   • RRPP package
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day5_gm_itv_annotated.R` (sections relevant to Q3) for context.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. SETUP — load packages and rebuild required objects
# -----------------------------------------------------------------------------
library(tidyverse)
library(geomorph)
library(RRPP)
library(lme4)
library(ggrepel)
library(patchwork)
library(FD)

theme_set(theme_bw(base_size = 12))
if (!dir.exists("outputs")) dir.create("outputs")

# Read landmark + metadata
coords <- readland.tps("data/fish_gm_day5.tps", specID = "ID", warnmsg = FALSE)
meta   <- read_csv("data/fish_population_metadata.csv", show_col_types = FALSE)
meta   <- meta[match(dimnames(coords)[[3]], meta$individual_id), ]

# GPA + shape PCA — required by all Day 5 questions
gpa <- gpagen(coords, print.progress = FALSE)
pca <- gm.prcomp(gpa$coords)

# =============================================================================
# Q3 — Population divergence and habitat effects (RRPP)
# =============================================================================
# -----------------------------------------------------------------------------
# Q3a — RRPP MANOVA: shape ~ habitat
# -----------------------------------------------------------------------------
# TASK
#   Fit `procD.lm(coords ~ habitat, data = gdf, iter = 999)`. Then fit
#   `procD.lm(coords ~ species + habitat, data = gdf)` to test the habitat
#   effect AFTER controlling for species. Interpret R² and p.
#
# EXPECTED OUTPUT
#   Two RRPP tables. Habitat effect partial (after species) is the
#   biologically meaningful one; expect R² of 0.05-0.15.
#
# HINTS
#   • Always include species as a covariate when testing habitat — otherwise
#     habitat absorbs species differences.
#   • Type-II SS via `anova(fit, type = 'II')` if available.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q3b — Pairwise Procrustes distances between populations
# -----------------------------------------------------------------------------
# TASK
#   Use `RRPP::pairwise()` to compute all pairwise distances between species ×
#   population groups. Plot as a heatmap. Are geographically close populations
#   (same river system) also morphologically close?
#
# EXPECTED OUTPUT
#   Heatmap outputs/group3_day5_Q3b_distances.png. Within-species pairs should
#   cluster; between-species distances dominate.
#
# HINTS
#   • `fit_pw <- pairwise(fit_pop, groups = paste(species, population,
#     sep='_'))`.
#   • `summary(fit_pw, test.type = 'dist')` returns the distance matrix.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q3c — Tangent space after removing species effect
# -----------------------------------------------------------------------------
# TASK
#   Extract residuals from the species-only model (`procD.lm(coords ~
#   species)`). Run a PCA on residuals. Plot — are populations separable on
#   the residual axes?
#
# EXPECTED OUTPUT
#   PNG outputs/group3_day5_Q3c_residuals_PCA.png. After removing species,
#   populations should partly cluster by habitat (fast vs slow).
#
# HINTS
#   • `fit_sp <- procD.lm(coords ~ species, data = gdf)`; `resid_arr <-
#     arrayspecs(fit_sp$residuals, p=12, k=2)`.
#   • PCA via `gm.prcomp(resid_arr)`.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group3_day5_*.png
# =============================================================================
