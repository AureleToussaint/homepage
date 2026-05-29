# =============================================================================
# DAY 5  ·  GROUP 2  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q2 exercise for GROUP 2.
# Submit your annotated R script with figures saved to outputs/ before Friday evening · final submission.
#
# RULES
#   • Use the ANNOTATED script (`day5_gm_itv_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group2_day5_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Quantify what fraction of body-shape variance lies between species, between populations within species, and within populations.
#
# DATA / OBJECTS YOU NEED
#   • gpa, pca_shape, meta (from Q1)
#   • lme4 package
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day5_gm_itv_annotated.R` (sections relevant to Q2) for context.
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
# Q2 — Variance partitioning of shape across hierarchical levels
# =============================================================================
# -----------------------------------------------------------------------------
# Q2a — Fit the hierarchical mixed model for PC1, PC2, PC3
# -----------------------------------------------------------------------------
# TASK
#   For each PC, fit `lmer(y ~ 1 + (1|species) + (1|species:population))`.
#   Extract variance components with `VarCorr()`.
#
# EXPECTED OUTPUT
#   Tibble with rows = PC and columns = grp (species / population / residual)
#   and pct (% variance).
#
# HINTS
#   • Use `as.data.frame(VarCorr(m))` and compute pct = vcov / sum(vcov).
#   • Loop with `lapply(1:3, function(i) ...)` and `bind_rows()`.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q2b — Stacked bar chart per PC
# -----------------------------------------------------------------------------
# TASK
#   Visualise the variance partition with a stacked bar (one bar per PC). Use
#   distinct colours for the three hierarchical levels.
#
# EXPECTED OUTPUT
#   PNG outputs/group2_day5_Q2b_variance.png. Between-species typically
#   dominates (60-75 %); between-population 15-25 %; residual 10-20 %.
#
# HINTS
#   • `geom_col(position = 'stack')`.
#   • Order factor levels: species, species:population, Residual.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q2c — Which axis has the highest between-population %?
# -----------------------------------------------------------------------------
# TASK
#   Identify the PC with the largest fraction of between-population variance
#   (relative to between-species). What biological process could explain this
#   — phenotypic plasticity, local adaptation, or ontogeny?
#
# EXPECTED OUTPUT
#   Identification + 5-line discussion in comments.
#
# HINTS
#   • Procrustes shape under fast-flow vs slow-flow water differs in
#     predictable ways: streamlined bodies, narrower peduncles → could load on
#     PC2 or PC3.
#   • Cite Riesch et al. (2022) for population-level shape divergence in
#     tropical fish.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group2_day5_*.png
# =============================================================================
