# =============================================================================
# DAY 3  ·  GROUP 3  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q3 exercise for GROUP 3.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 4.
#
# RULES
#   • Use the ANNOTATED script (`day3_rlq_fourthcorner_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group3_day3_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Add a 4th extinction scenario (most-common first), test whether removing only the CR + EN species causes more functional erosion than expected at random, and locate the species with the highest functional uniqueness.
#
# DATA / OBJECTS YOU NEED
#   • outputs/day3_extinction_curves.rds
#   • outputs/day3_objects.rds (trait_mat_sc, comm_mat, fuse_df)
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day3_rlq_fourthcorner_annotated.R` (sections relevant to Q3) for context.
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
# Q3 — Extinction simulations: scenarios and quantitative tests
# =============================================================================
# -----------------------------------------------------------------------------
# Q3a — Add 4th scenario: most common first
# -----------------------------------------------------------------------------
# TASK
#   Build a sequence ordering species from most common to least common
#   (inverse of scenario B). Run the same `fd_along()` loop. Plot all 4
#   scenarios on the same figure.
#
# EXPECTED OUTPUT
#   PNG outputs/group3_day3_Q3a_4scenarios.png. Discuss: removing common
#   species first usually preserves FRic longer than removing rare ones, IF
#   the rare species are functionally unique.
#
# HINTS
#   • `seq_common <- names(sort(colSums(comm_mat), decreasing = TRUE))`.
#   • Reuse the `fd_along()` function from the lecture script.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q3b — Empirical test: %FRic loss with CR + EN only
# -----------------------------------------------------------------------------
# TASK
#   Identify CR and EN species from `iucn`. Compute the observed %FRic lost
#   when only these are removed. Compare against 199–999 random removals of
#   the same number of species. Compute an empirical p-value (one-sided).
#
# EXPECTED OUTPUT
#   Numeric output: observed loss, mean random loss, empirical p. Density plot
#   showing the random distribution and the observed value as a vertical line.
#
# HINTS
#   • Wrap the random sample in `replicate(199, ...)` for speed; bump to 999
#     for the final analysis.
#   • Empirical p = `mean(loss_perm >= observed)`.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q3c — Top-FUn species in PCoA space
# -----------------------------------------------------------------------------
# TASK
#   Find the species with the highest normalised FUn from `fuse_df`. Plot the
#   PCoA of all species; highlight the top-FUn species and label its 2 nearest
#   neighbours. What makes it unique?
#
# EXPECTED OUTPUT
#   PNG outputs/group3_day3_Q3c_topFUn.png. The species (often Arapaima or
#   Electrophorus) sits at the periphery of the trait cloud.
#
# HINTS
#   • `pcoa <- ape::pcoa(dist(trait_mat_sc))`.
#   • Use a different colour or larger point size for the focal species.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group3_day3_*.png
# =============================================================================
