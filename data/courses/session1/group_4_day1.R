# =============================================================================
# DAY 1  ·  GROUP 4  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q4 exercise for GROUP 4.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 2.
#
# RULES
#   • Use the ANNOTATED script (`day1_trait_data_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group4_day1_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Quantify within-species variability of SLA using TRY-format observations; compare ITV across PFTs; critically discuss limitations of TRY-based ITV estimates.
#
# DATA / OBJECTS YOU NEED
#   • try_filtered  (long-format trait table loaded by Day 1 script, with columns species, trait_short, StdValue, PFT_label)
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day1_trait_data_annotated.R` (sections relevant to Q4) for context.
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
# Q4 — Intraspecific trait variability (ITV) from TRY long-format
# =============================================================================
# -----------------------------------------------------------------------------
# Q4a — CV of SLA per species
# -----------------------------------------------------------------------------
# TASK
#   Filter try_filtered to SLA observations only. Compute the coefficient of
#   variation (CV = sd / mean) per species. Identify the species with the
#   highest CV.
#
# EXPECTED OUTPUT
#   Tibble with one row per species and columns mean, sd, CV. The top species
#   should be flagged with print() or sliced via `slice_max`.
#
# HINTS
#   • `filter(trait_short == 'SLA') |> group_by(species)`.
#   • Use `summarise(mean = mean(StdValue, na.rm = TRUE), sd = sd(StdValue,
#     na.rm = TRUE))` then mutate CV.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q4b — CV distribution by PFT
# -----------------------------------------------------------------------------
# TASK
#   Boxplot of CVs grouped by PFT. Which PFT shows the highest median CV? Test
#   the difference across PFTs (Kruskal-Wallis).
#
# EXPECTED OUTPUT
#   PNG outputs/group4_day1_Q4b_CV_PFT.png + Kruskal-Wallis p-value. Ruderals
#   often show high CV (large environmental responsiveness).
#
# HINTS
#   • `kruskal.test(CV ~ PFT_label, data = ...)` — non-parametric ANOVA.
#   • Beware low n per PFT — discuss this when interpreting the p-value.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q4c — Correlation between mean SLA and CV
# -----------------------------------------------------------------------------
# TASK
#   Compute Spearman correlation between species mean SLA and species CV. Plot
#   the relationship. Do high-SLA species show proportionally larger or
#   smaller variability?
#
# EXPECTED OUTPUT
#   A Spearman r value, p-value and a labelled scatter plot. A negative
#   correlation would suggest stress-tolerators are more variable; positive
#   would suggest acquisitive species are more variable.
#
# HINTS
#   • `cor(mean_SLA, CV, method = 'spearman')`.
#   • `geom_smooth(method = 'lm')` for visual aid.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q4d — Critical discussion of TRY-based ITV
# -----------------------------------------------------------------------------
# TASK
#   Write a 6-10 line critique answering: does the CV computed here really
#   capture biological ITV? What confounds it (sites, age, measurement
#   protocols)? How would a common-garden experiment differ?
#
# EXPECTED OUTPUT
#   Comment block at the end of the script.
#
# HINTS
#   • TRY pools measurements from many sites/years/observers — confounds ITV
#     (genuine within-species variation) with measurement noise and between-
#     population differences.
#   • Common-garden experiments isolate the genetic component by removing the
#     environmental gradient.
#   • Cite Violle et al. (2012) TREE on the importance of ITV partitioning.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group4_day1_*.png
# =============================================================================
