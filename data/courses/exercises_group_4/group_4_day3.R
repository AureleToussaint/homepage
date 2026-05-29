# =============================================================================
# DAY 3  ·  GROUP 4  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q4 exercise for GROUP 4.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 4.
#
# RULES
#   • Use the ANNOTATED script (`day3_rlq_fourthcorner_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group4_day3_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Identify the top-5 FUSE conservation priorities for Amazonian fish, design an Amazon-specific re-weighting of IUCN categories, and discuss WHY functionally unique species are more vulnerable.
#
# DATA / OBJECTS YOU NEED
#   • outputs/day3_fuse_scores.rds
#   • outputs/day3_extinction_curves.rds
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day3_rlq_fourthcorner_annotated.R` (sections relevant to Q4) for context.
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
# Q4 — FUSE for Amazonian fish: priorities and regional adaptation
# =============================================================================
# -----------------------------------------------------------------------------
# Q4a — Top-5 conservation priorities
# -----------------------------------------------------------------------------
# TASK
#   From `fuse_df`, extract the top-5 species by FUSE score. Compare with the
#   first 5 species lost under scenario C (most threatened first) — how much
#   overlap?
#
# EXPECTED OUTPUT
#   CSV outputs/group4_day3_Q4a_top5.csv with species, IUCN, FUn, FSpec and
#   FUSE. Overlap count printed.
#
# HINTS
#   • `slice_max(fuse_df, FUSE, n = 5)`.
#   • Build the threat-sorted sequence from `iucn_w` and intersect.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q4b — Amazon-specific IUCN weights
# -----------------------------------------------------------------------------
# TASK
#   Propose alternative IUCN weights tailored to Amazonian deforestation
#   pressures. Justify each weight in 1-2 lines. Recompute FUSE with these
#   weights and compare the top-5 lists.
#
# EXPECTED OUTPUT
#   CSV outputs/group4_day3_Q4b_amazon_weights.csv with original and Amazon
#   weights and a one-line rationale per category.
#
# HINTS
#   • Suggested: bump LC weight to 0.05 (deforestation pressure on common
#     spp), VU to 0.50, DD to 0.50 (precautionary).
#   • Recompute FUSE and compare rankings with Spearman correlation.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q4c — Discussion ≥ 150 words
# -----------------------------------------------------------------------------
# TASK
#   Why are functionally unique species more vulnerable than redundant ones?
#   Address: (i) ecological mechanisms (small populations, narrow niches, low
#   fecundity); (ii) life-history correlates (large body, late maturity);
#   (iii) implications for ecosystem functioning if these species are lost.
#   Cite ≥ 2 references.
#
# EXPECTED OUTPUT
#   Comment block at the end of the script.
#
# HINTS
#   • Mouillot et al. (2013) PLoS Biol — rare species support vulnerable
#     functions.
#   • Carmona et al. (2021) Sci Adv — global functional erosion.
#   • Slow-pace life-history → lower demographic resilience.
#   • Functional uniqueness without redundancy = no insurance for functions.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group4_day3_*.png
# =============================================================================
