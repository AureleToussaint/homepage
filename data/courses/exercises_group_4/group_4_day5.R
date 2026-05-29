# =============================================================================
# DAY 5  ·  GROUP 4  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q4 exercise for GROUP 4.
# Submit your annotated R script with figures saved to outputs/ before Friday evening · final submission.
#
# RULES
#   • Use the ANNOTATED script (`day5_gm_itv_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group4_day5_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Compare community-level FRic computed from the 10 TM trait ratios (Day 4) and from the GM shape PC scores (Day 5). Discuss when ITV must be integrated into FD analyses.
#
# DATA / OBJECTS YOU NEED
#   • data/fish_traits_sp_means.csv (TM)
#   • data/fish_communities_amazon.csv
#   • gpa, pca, meta from Q1
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day5_gm_itv_annotated.R` (sections relevant to Q4) for context.
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
# Q4 — Course synthesis: TM-based vs GM-based functional rankings
# =============================================================================
# -----------------------------------------------------------------------------
# Q4a — FRic from TM traits (10 ratios)
# -----------------------------------------------------------------------------
# TASK
#   Compute FRic per Amazonian site using the scaled 10-trait matrix (Day 4
#   dataset).
#
# EXPECTED OUTPUT
#   Tibble of 12 sites × FRic_TM. The same as Day 4 Q3a — re-do for comparison
#   with GM.
#
# HINTS
#   • Restrict to species present in BOTH the TM and GM datasets.
#   • Use `dbFD(..., m = 4, stand.x = FALSE)`.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q4b — FRic from GM shape PCs
# -----------------------------------------------------------------------------
# TASK
#   Aggregate the 240 specimens to 5 species means in shape PC space (PC1-3).
#   Compute FRic per site using these 3 axes (only the 5 GM species present in
#   the community matrix).
#
# EXPECTED OUTPUT
#   Tibble of 12 sites × FRic_GM (computed only on the 5 GM species). Compare
#   TM and GM with Spearman rank correlation.
#
# HINTS
#   • `sp_means_pc <- pca$x[,1:3] |> as.data.frame() |>
#     rownames_to_column('individual_id') |> left_join(meta) |>
#     group_by(species) |> summarise(across(starts_with('Comp'), mean))`.
#   • Restrict community matrix to the 5 GM species.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q4c — When does ignoring ITV bias FD estimates?
# -----------------------------------------------------------------------------
# TASK
#   Write a 150+ word discussion answering: WHEN does ignoring intraspecific
#   variability most strongly bias FD? Use Albert et al. (2010) and Violle et
#   al. (2012) as anchors. Apply to Amazonian fish: should we recommend ITV-
#   aware FD by default?
#
# EXPECTED OUTPUT
#   Comment block at the end of the script.
#
# HINTS
#   • Bias is largest when (i) within-species variance is large relative to
#     between-species, (ii) populations span the environmental gradient
#     asymmetrically, (iii) species are filtered along that same gradient.
#   • Albert et al. (2010) Funct. Ecol. — partition of variance test.
#   • Practically: ITV matters most for between-population comparisons and for
#     low-richness communities.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group4_day5_*.png
# =============================================================================
