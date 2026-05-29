# =============================================================================
# DAY 3  ·  GROUP 1  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q1 exercise for GROUP 1.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 4.
#
# RULES
#   • Use the ANNOTATED script (`day3_rlq_fourthcorner_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group1_day3_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Interpret RLQ outputs robustly: test how a single trait removal changes the main axis, identify driver species, and check for spatial autocorrelation in the trait–environment relationship.
#
# DATA / OBJECTS YOU NEED
#   • outputs/day3_objects.rds (R, L, Q, trait_mat_sc, comm_mat, fuse_df)
#   • outputs/day3_rlq_result.rds
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day3_rlq_fourthcorner_annotated.R` (sections relevant to Q1) for context.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. SETUP — load packages and rebuild required objects
# -----------------------------------------------------------------------------
library(tidyverse)
library(ade4)
library(FD)
library(funrar)
library(ggrepel)
library(patchwork)
library(ape)

theme_set(theme_bw(base_size = 12))
if (!dir.exists("outputs")) dir.create("outputs")

# Run the annotated Day 3 script first; it saves the objects we need below.
o            <- readRDS("outputs/day3_objects.rds")
fish_traits  <- read_csv("data/fish_traits_sp_means.csv",     show_col_types = FALSE)
fish_comm    <- read_csv("data/fish_communities_amazon.csv",  show_col_types = FALSE)
site_env     <- read_csv("data/fish_site_env.csv",            show_col_types = FALSE)
iucn         <- read_csv("data/fish_iucn_status.csv",         show_col_types = FALSE)
rlq_res      <- readRDS("outputs/day3_rlq_result.rds")
fc_combined  <- readRDS("outputs/day3_fourthcorner_combined.rds")
fc_tab       <- readRDS("outputs/day3_fc_heatmap_data.rds")
ext_all      <- readRDS("outputs/day3_extinction_curves.rds")
fuse_df      <- readRDS("outputs/day3_fuse_scores.rds")

R_tab <- o$R; L_tab <- o$L; Q_tab <- o$Q
trait_mat_sc <- o$trait_mat_sc; comm_mat <- o$comm_mat

# =============================================================================
# Q1 — RLQ interpretation and spatial autocorrelation
# =============================================================================
# -----------------------------------------------------------------------------
# Q1a — Re-run RLQ excluding Body_size from Q
# -----------------------------------------------------------------------------
# TASK
#   Build `Q_noSize` (Q without Body_size). Re-run the chain dudi.coa →
#   weighted PCAs → rlq(). Compare Axis 1 environmental loadings to those of
#   the full model. Plot side-by-side bars.
#
# EXPECTED OUTPUT
#   Bar plot outputs/group1_day3_Q1a_rlq_noBodySize.png and a 3-line
#   interpretation: does removing Body_size change which environmental
#   gradient dominates Axis 1?
#
# HINTS
#   • Weighted PCAs need `row.w = dudiL$lw` for R and `row.w = dudiL$cw` for
#     Q.
#   • Compare `rlq_full$l1[,1]` and `rlq_noSize$l1[,1]`.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q1b — Top-5 species at each end of Axis 1
# -----------------------------------------------------------------------------
# TASK
#   Sort species by their RLQ Axis 1 score. Extract the 5 highest and 5
#   lowest. Bar-plot them. List which traits dominate at each extreme by
#   linking back to the trait loadings c1[,1].
#
# EXPECTED OUTPUT
#   Bar plot outputs/group1_day3_Q1b_extreme_species.png. Comment on which
#   trait combinations characterise each extreme.
#
# HINTS
#   • `rlq_res$lQ` contains species scores.
#   • Use `slice_max` / `slice_min` after converting to a tibble.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q1c — Moran's I on CWM(Body_elong) ~ degradation residuals
# -----------------------------------------------------------------------------
# TASK
#   Compute CWM(Body_elong) per site. Fit `lm(CWM_BE ~ degradation)`. Compute
#   Moran's I on the residuals using site lat/lon coordinates. Is there
#   spatial structure remaining after the environmental covariate is accounted
#   for?
#
# EXPECTED OUTPUT
#   Moran's I value, p-value, and a brief interpretation. p < 0.05 suggests
#   spatial autocorrelation that should be modelled (e.g. add a spatial random
#   effect).
#
# HINTS
#   • Build inverse-distance weights: `W <- 1 / as.matrix(dist(coords))`; set
#     diagonal to 0; row-normalise.
#   • `ape::Moran.I(residuals, W)` returns observed I, expected and p-value.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group1_day3_*.png
# =============================================================================
