# =============================================================================
# DAY 4  ·  GROUP 4  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q4 exercise for GROUP 4.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 5.
#
# RULES
#   • Use the ANNOTATED script (`day4_fish_landmarks_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group4_day4_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Compute FUSE scores for the 30 Amazonian species using the 10 morphological traits, identify the most irreplaceable species, and contrast the morphology-based fish approach with the leaf-based plant approach of Days 1-3.
#
# DATA / OBJECTS YOU NEED
#   • data/fish_traits_sp_means.csv
#   • data/fish_iucn_status.csv
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day4_fish_landmarks_annotated.R` (sections relevant to Q4) for context.
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
# Q4 — From morphology to function: FUSE for Amazonian fish
# =============================================================================
# -----------------------------------------------------------------------------
# Q4a — Compute FUSE for the 30 fish species
# -----------------------------------------------------------------------------
# TASK
#   Implement FUSE = FUn_norm × FSpec_norm × IUCN_weight using the scaled
#   10-trait matrix. Use Mooers et al. (2008) global weights.
#
# EXPECTED OUTPUT
#   Tibble sorted by FUSE descending; top species printed.
#
# HINTS
#   • FUn = nearest-neighbour distance in trait space.
#   • FSpec = distance from centroid.
#   • Normalise both to [0, 1] before multiplying by IUCN weight.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q4b — Bubble chart of FUSE
# -----------------------------------------------------------------------------
# TASK
#   Plot FUn_norm (x) × FSpec_norm (y), bubble size proportional to FUSE,
#   colour by IUCN status. Label the top-3 species.
#
# EXPECTED OUTPUT
#   PNG outputs/group4_day4_Q4b_FUSE.png. Top species (Arapaima,
#   Pseudoplatystoma, Brachyplatystoma) sit in the upper-right corner.
#
# HINTS
#   • `scale_size(range = c(2, 9))`.
#   • `ggrepel::geom_text_repel(data = top3, ...)` for selective labelling.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q4c — Discussion: fish morphology vs plant leaves (≥ 150 words)
# -----------------------------------------------------------------------------
# TASK
#   Contrast the trait-based approach for fish (morphological traits from
#   photographs) and plants (leaf economics traits from field). Address: ease
#   of measurement, link to function, taxonomic universality, limitations
#   specific to each system.
#
# EXPECTED OUTPUT
#   Comment block at end of script.
#
# HINTS
#   • Plant traits (SLA, LDMC) directly measure C/water economy; fish traits
#     indirectly proxy locomotion and diet.
#   • Plant traits require fresh material; fish traits can be measured from
#     preserved specimens or photographs.
#   • Both standardised globally (TRY for plants, FishMorph for fish).

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group4_day4_*.png
# =============================================================================
