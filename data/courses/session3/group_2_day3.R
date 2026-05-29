# =============================================================================
# DAY 3  ·  GROUP 2  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q2 exercise for GROUP 2.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 4.
#
# RULES
#   • Use the ANNOTATED script (`day3_rlq_fourthcorner_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group2_day3_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Identify the strongest trait × environment associations, propose mechanistic hypotheses, and demonstrate the Type-I inflation that occurs when only one permutation model is used.
#
# DATA / OBJECTS YOU NEED
#   • outputs/day3_fc_heatmap_data.rds (combined max-test)
#   • R_tab, L_tab, Q_tab from day3_objects.rds
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day3_rlq_fourthcorner_annotated.R` (sections relevant to Q2) for context.
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
# Q2 — Critical assessment of the fourth-corner test
# =============================================================================
# -----------------------------------------------------------------------------
# Q2a — Three strongest significant associations
# -----------------------------------------------------------------------------
# TASK
#   From `fc_tab` (loaded), filter to p.adj < 0.05 and slice the 3 rows with
#   largest |stat|. For EACH of the three, propose a mechanistic biological
#   hypothesis citing one published study.
#
# EXPECTED OUTPUT
#   Tibble of 3 rows + 3 paragraphs (4-6 lines each) of biological hypotheses
#   with citations.
#
# HINTS
#   • Body_elong × current_velocity → hydrodynamic drag (Villéger et al. 2010
#     Ecol. Appl.).
#   • Pec_pos × turbidity → benthic specialisation (Toussaint et al. 2016 Sci
#     Rep).
#   • Body_size × degradation → fishing pressure on large-bodied species (Su
#     et al. 2021 Science).

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q2b — Type-I inflation: model 2 alone vs max-test
# -----------------------------------------------------------------------------
# TASK
#   Re-run `fourthcorner()` with `modeltype = 2` ONLY (no max-test). Count
#   significant cells (p.adj < 0.05) in this version vs the combined max-test.
#   Bar plot the two counts.
#
# EXPECTED OUTPUT
#   Bar plot outputs/group2_day3_Q2b_inflation.png. Model 2 alone produces ~2×
#   more significant cells — illustrates Type-I inflation.
#
# HINTS
#   • Don't forget `set.seed(42)` for reproducibility.
#   • `fc_m2$tabD2$adj.pvalue < 0.05`.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q2c — Linearity check on the strongest pair
# -----------------------------------------------------------------------------
# TASK
#   Take the strongest significant pair (Q2a). Plot the corresponding site-
#   level scatter (CWM(trait) ~ environment). Fit linear and quadratic LMs.
#   Compare with `anova()` likelihood-ratio. Is the linear assumption
#   justified?
#
# EXPECTED OUTPUT
#   Scatter with two superimposed smoothers (linear + quadratic), anova
#   p-value. Discuss whether non-linearity matters for the fourth-corner
#   conclusions.
#
# HINTS
#   • `anova(mod_lin, mod_quad)` gives F and p for the additional term.
#   • Even significant non-linearity may not invalidate the fourth-corner
#     result; it just refines the response shape.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group2_day3_*.png
# =============================================================================
