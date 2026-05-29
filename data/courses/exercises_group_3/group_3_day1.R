# =============================================================================
# DAY 1  ·  GROUP 3  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q3 exercise for GROUP 3.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 2.
#
# RULES
#   • Use the ANNOTATED script (`day1_trait_data_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group3_day1_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Locate the most acquisitive and most conservative species of the dataset in the LES space; build a multi-trait composite score to rank species by their position along the fast–slow continuum.
#
# DATA / OBJECTS YOU NEED
#   • sp_traits  (with SLA, LDMC, LeafN, Height, SeedMass, PFT_label, IUCN_status)
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day1_trait_data_annotated.R` (sections relevant to Q3) for context.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. SETUP — load packages and rebuild required objects
# -----------------------------------------------------------------------------
library(tidyverse)
library(ggrepel)
library(patchwork)
library(GGally)

theme_set(theme_bw(base_size = 12))
if (!dir.exists("outputs")) dir.create("outputs")

# Run the annotated script first to load sp_traits, try_filtered, ...
source("day1_trait_data_annotated.R")

# =============================================================================
# Q3 — Extreme species in the LES space and a composite acquisitiveness score
# =============================================================================
# -----------------------------------------------------------------------------
# Q3a — Identify SLA extremes
# -----------------------------------------------------------------------------
# TASK
#   Identify the species with minimum and maximum SLA, plus the top-3 and
#   bottom-3. Report their PFT, IUCN status and other trait values.
#
# EXPECTED OUTPUT
#   Two short tables (top-3 highest SLA, bottom-3 lowest SLA) with the full
#   trait profile per species.
#
# HINTS
#   • `slice_max(SLA, n = 3)` and `slice_min(SLA, n = 3)`.
#   • Expected: SLA min in PFT_stress (Dryas, Carex_firma); SLA max in
#     PFT_ruderal (Rumex, Chenopodium).

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q3b — Position in LES bivariate plots
# -----------------------------------------------------------------------------
# TASK
#   Produce two LES biplots: log(SLA) × log(LeafN) and log(SLA) × log(LDMC).
#   Highlight the SLA extremes (different colour or shape). Comment on their
#   position relative to the LES diagonal.
#
# EXPECTED OUTPUT
#   Two-panel figure saved as outputs/group3_day1_Q3b_LES.png. The low-SLA
#   species cluster in the conservative corner; high-SLA in the acquisitive
#   corner.
#
# HINTS
#   • Use `mutate(highlight = species %in% c(...))` to flag extremes.
#   • `ggrepel::geom_text_repel` with `data = filter(., highlight)`.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q3c — Multi-trait coherence of the SLA-max species
# -----------------------------------------------------------------------------
# TASK
#   Test whether the species with maximum SLA also has (i) lowest SeedMass,
#   (ii) highest Height. Discuss WHY these expectations hold or fail (traits
#   do NOT all line up perfectly).
#
# EXPECTED OUTPUT
#   Comparison table; explanatory comment about LES vs other functional axes
#   (Díaz et al. 2016 found ≥ 2 orthogonal axes globally).
#
# HINTS
#   • Pull the SLA-max row; check rank of its other traits across all spp.
#   • Diaz et al. 2016 Nature: size axis (Height, SeedMass) is orthogonal to
#     LES (SLA, LeafN) — perfect alignment is not expected.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q3d — Composite acquisitiveness score
# -----------------------------------------------------------------------------
# TASK
#   For each species, compute the rank (1 = lowest) of its SLA, LeafN and
#   inverted LDMC values. Sum the three to obtain a composite score. Identify
#   the most 'acquisitive' species (highest score). Plot rank-rank panels.
#
# EXPECTED OUTPUT
#   Tibble of species sorted by composite score. The top species should belong
#   to PFT_ruderal.
#
# HINTS
#   • `mutate(rank_SLA = rank(SLA), rank_LeafN = rank(LeafN), rank_LDMC_inv =
#     rank(-LDMC))`.
#   • Composite = sum of the three ranks. Higher = more acquisitive.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group3_day1_*.png
# =============================================================================
