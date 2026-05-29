# =============================================================================
# DAY 4  ·  GROUP 2  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q2 exercise for GROUP 2.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 5.
#
# RULES
#   • Use the ANNOTATED script (`day4_fish_landmarks_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group2_day4_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Visualise how the 10 morphological traits separate trophic guilds, identify trait-pair correlations and locate the central / peripheral species in PCoA space.
#
# DATA / OBJECTS YOU NEED
#   • data/fish_traits_sp_means.csv
#   • data/fish_species_metadata.csv (trophic_guild, family, order)
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day4_fish_landmarks_annotated.R` (sections relevant to Q2) for context.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. SETUP — load packages and rebuild required objects
# -----------------------------------------------------------------------------
library(tidyverse)
library(FD)
library(vegan)
library(ggrepel)
library(ggforce)
library(patchwork)
library(corrplot)
library(irr)

theme_set(theme_bw(base_size = 12))
if (!dir.exists("outputs")) dir.create("outputs")

# All Day 4 objects load directly from data/ — no Day 3 dependency.

# =============================================================================
# Q2 — Trait–ecology relationships and the Neotropical functional space
# =============================================================================
# -----------------------------------------------------------------------------
# Q2a — Mouth_pos × Body_elong by trophic guild
# -----------------------------------------------------------------------------
# TASK
#   Scatter Mouth_pos vs Body_elong, colour by trophic guild. Are guilds
#   visually separable? Add ellipses or convex hulls per guild.
#
# EXPECTED OUTPUT
#   PNG outputs/group2_day4_Q2a_trait_space.png. Surface predators should
#   occupy high Mouth_pos × high Body_elong; benthic feeders low Mouth_pos.
#
# HINTS
#   • `stat_ellipse(aes(group = trophic_guild), level = 0.68)`.
#   • Or `ggforce::geom_mark_hull()` for convex hulls.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q2b — Full trait correlation matrix
# -----------------------------------------------------------------------------
# TASK
#   Compute the 10 × 10 Pearson correlation matrix of all traits. Visualise
#   with `corrplot()`. Identify the most strongly correlated trait pair and
#   propose a mechanistic explanation (locomotion–diet trade-off).
#
# EXPECTED OUTPUT
#   PNG outputs/group2_day4_Q2b_trait_corr.png + 1-paragraph interpretation of
#   the strongest pair.
#
# HINTS
#   • Body_size × Oral_gape often r > 0.8 — large fish have large gapes
#     (allometric constraint).
#   • Body_elong × Caud_AR often correlated — both reflect hydrodynamic
#     adaptation.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q2c — PCoA of the 10-trait functional space
# -----------------------------------------------------------------------------
# TASK
#   Compute Euclidean distances on scaled traits, run `ape::pcoa()`. Plot Axes
#   1-2, colour by trophic guild OR by family. Which order of fish occupies
#   the most extreme position? Which is the most central?
#
# EXPECTED OUTPUT
#   PNG outputs/group2_day4_Q2c_PCoA.png. Siluriformes (benthic) and
#   Osteoglossiformes (Arapaima) occupy peripheral positions; Characiformes is
#   the most central order.
#
# HINTS
#   • Use `ape::pcoa(dist(scale(trait_mat)))`.
#   • Plot variance explained by PC1 / PC2 in the axis labels.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group2_day4_*.png
# =============================================================================
