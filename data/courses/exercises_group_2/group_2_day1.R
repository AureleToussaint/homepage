# =============================================================================
# DAY 1  ·  GROUP 2  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q2 exercise for GROUP 2.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 2.
#
# RULES
#   • Use the ANNOTATED script (`day1_trait_data_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group2_day1_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Quantify and interpret the allometric relationship between seed mass and plant height; demonstrate why log-transformation is statistically and biologically necessary; test for PFT-specific allometric scaling.
#
# DATA / OBJECTS YOU NEED
#   • sp_traits  (with SeedMass, Height, PFT_label)
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day1_trait_data_annotated.R` (sections relevant to Q2) for context.
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
# Q2 — Allometry: SeedMass × Height and the role of log-transformation
# =============================================================================
# -----------------------------------------------------------------------------
# Q2a — Pearson correlation: raw vs log scales
# -----------------------------------------------------------------------------
# TASK
#   Compute Pearson correlation between SeedMass and Height on (i) raw values,
#   (ii) log-transformed values. Test the significance of the log-log
#   correlation with `cor.test()`.
#
# EXPECTED OUTPUT
#   Two correlation values (typically r_log > r_raw), and a cor.test p-value.
#   The log-log relationship should be much stronger.
#
# HINTS
#   • `cor(x, y, use = 'complete.obs')` for the coefficient.
#   • `cor.test(log(x), log(y))` for the inferential test.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q2b — Two-panel scatterplot (raw + log-log)
# -----------------------------------------------------------------------------
# TASK
#   Produce a side-by-side comparison: left panel = raw scale, right panel =
#   log-log. Add the regression line and the correlation value in the
#   subtitle. Label the species on the log-log panel.
#
# EXPECTED OUTPUT
#   PNG saved as outputs/group1_day1_Q2b_allometry.png. The log-log panel
#   should show a clear linear cloud.
#
# HINTS
#   • Use `ggrepel::geom_text_repel()` for non-overlapping species labels.
#   • Combine panels with `patchwork`.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q2c — Justification of log-transformation
# -----------------------------------------------------------------------------
# TASK
#   Plot histograms of the four variables (raw and log) and write a comment
#   block (statistical AND biological justification) explaining why we work in
#   log space.
#
# EXPECTED OUTPUT
#   2×2 histogram panel; ~10 lines of justification text in comments.
#
# HINTS
#   • Statistical: log normalises right-skewed distributions, equalises weight
#     of small and large species.
#   • Biological: allometric scaling laws (West et al. 1997 Science) follow
#     power laws → linear in log-log space; the slope is biologically
#     interpretable.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q2d — Heterogeneous slopes between PFTs?
# -----------------------------------------------------------------------------
# TASK
#   Fit two models: (i) `lm(log(Height) ~ log(SeedMass) * PFT_label)` with
#   interaction; (ii) without interaction. Compare with AIC. Test the global
#   interaction with `anova()`.
#
# EXPECTED OUTPUT
#   ANOVA table, AIC comparison, and a faceted plot showing per-PFT regression
#   lines. Discuss the result given the small n per PFT.
#
# HINTS
#   • With n ≈ 6 species per PFT, statistical power is limited — discuss this
#     caveat regardless of the p-value.
#   • Visual inspection of slopes is as informative as the formal test.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group2_day1_*.png
# =============================================================================
