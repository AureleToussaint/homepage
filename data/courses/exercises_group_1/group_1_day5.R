# =============================================================================
# DAY 5  ·  GROUP 1  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q1 exercise for GROUP 1.
# Submit your annotated R script with figures saved to outputs/ before Friday evening · final submission.
#
# RULES
#   • Use the ANNOTATED script (`day5_gm_itv_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group1_day5_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Run a Generalised Procrustes Analysis on a multi-population fish dataset, interpret the shape PCs through deformation grids, and test for allometry.
#
# DATA / OBJECTS YOU NEED
#   • data/fish_gm_day5.tps
#   • data/fish_population_metadata.csv
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day5_gm_itv_annotated.R` (sections relevant to Q1) for context.
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
# Q1 — Procrustes PCA, deformation grids and allometry
# =============================================================================
# -----------------------------------------------------------------------------
# Q1a — Run GPA and shape PCA
# -----------------------------------------------------------------------------
# TASK
#   Read the .tps file with `readland.tps()`. Run `gpagen()`. Run
#   `gm.prcomp()`. Plot PC1 × PC2 coloured by species, shaped by population.
#   Compute % variance for the first 5 PCs.
#
# EXPECTED OUTPUT
#   PNG outputs/group1_day5_Q1a_shape_PCA.png with cluster ellipses by
#   population × species. PC1 typically captures 30-50 % of shape variation.
#
# HINTS
#   • Verify dim(coords) is p=12, k=2, n=240.
#   • Use `stat_ellipse(aes(group = paste(species, population)))`.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q1b — Deformation grids along PC1
# -----------------------------------------------------------------------------
# TASK
#   Compute the consensus shape with `mshape()`. Use `plotRefToTarget()` to
#   show the shape deformation between the consensus and the minimum / maximum
#   of PC1 (thin-plate spline).
#
# EXPECTED OUTPUT
#   PNG outputs/group1_day5_Q1b_deformation.png with two side-by-side TPS
#   grids.
#
# HINTS
#   • Use `pca_shape$shapes$shapes.comp1$min` and `$max`.
#   • Set `mag = 2` to amplify visualisation.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q1c — Allometry test
# -----------------------------------------------------------------------------
# TASK
#   Test whether log(centroid size) explains shape: `procD.lm(coords ~
#   log(Csize), data = gdf, iter = 999)`. Plot PC1 × log(Csize) coloured by
#   species. Discuss: does body size predict shape across the whole dataset?
#
# EXPECTED OUTPUT
#   RRPP-MANOVA p-value, R² + scatter plot. R² typically 0.10-0.30 across
#   species; can be much higher within species.
#
# HINTS
#   • Use `geomorph.data.frame()` to build the gdf.
#   • Within-species allometry: stratify by species using the same `procD.lm`
#     with a `~ log(Csize) * species` interaction.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group1_day5_*.png
# =============================================================================
