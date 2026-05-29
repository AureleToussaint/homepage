# =============================================================================
# DAY 1  ·  GROUP 1  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q1 exercise for GROUP 1.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 2.
#
# RULES
#   • Use the ANNOTATED script (`day1_trait_data_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group1_day1_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Test whether IUCN-threatened species in the alpine grassland dataset have different SLA values from non-threatened species, and interpret the result in light of leaf economics theory.
#
# DATA / OBJECTS YOU NEED
#   • sp_traits  (loaded by day1_trait_data_annotated.R)
#   •   Required columns: species, SLA, IUCN_status, threatened, PFT_label
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day1_trait_data_annotated.R` (sections relevant to Q1) for context.
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
# Q1 — IUCN status and SLA: are threatened species more conservative?
# =============================================================================
# -----------------------------------------------------------------------------
# Q1a — Compute the % of threatened species
# -----------------------------------------------------------------------------
# TASK
#   Calculate the proportion of species classified as VU, EN or CR (the
#   logical column `threatened` is already coded). Also produce a count table
#   of the 5 IUCN categories with their relative percentage.
#
# EXPECTED OUTPUT
#   A single value (≈ 35–45 %) and a small table with one row per IUCN
#   category ordered LC → CR.
#
# HINTS
#   • Use `summarise()` with `n()` and `sum(threatened)`.
#   • For the count table: `count(IUCN_status)` then mutate a percentage
#     column.
#   • Re-order the factor levels with `factor(..., levels =
#     c('LC','NT','VU','EN','CR'))`.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q1b — Visualise SLA distribution: threatened vs non-threatened
# -----------------------------------------------------------------------------
# TASK
#   Produce a violin + boxplot comparing SLA between the two groups, colouring
#   points by PFT. Optionally add a density plot on log(SLA).
#
# EXPECTED OUTPUT
#   Two-panel figure (violin/box on left, density on right) saved as
#   outputs/group1_day1_Q1b_SLA_IUCN.png. Threatened species typically shift
#   toward LOWER SLA values.
#
# HINTS
#   • Use `geom_violin()` + `geom_boxplot(width = 0.15)` + `geom_jitter()`.
#   • Map `fill = threatened` and `colour = PFT_label`.
#   • Combine the two plots with `patchwork::+` operator.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q1c — Statistical test (Wilcoxon)
# -----------------------------------------------------------------------------
# TASK
#   Test whether the SLA distribution differs between threatened and non-
#   threatened species. Justify why a Wilcoxon test (non-parametric) is
#   preferable to a Student t-test here.
#
# EXPECTED OUTPUT
#   A `wilcox.test` output with W statistic and p-value. Report median, mean,
#   sd and IQR per group.
#
# HINTS
#   • `wilcox.test(SLA ~ threatened, data = sp_traits, exact = FALSE, conf.int
#     = TRUE)`
#   • Justification: small n per group (~ 12–18), no normality assumption
#     needed.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q1d — Ecological interpretation
# -----------------------------------------------------------------------------
# TASK
#   Write a 5-7 line interpretation answering: (i) is SLA higher or lower in
#   threatened species? (ii) which leaf-economics strategy does this
#   correspond to? (iii) which PFTs drive this pattern? (iv) cite one
#   supporting reference from the course bibliography.
#
# EXPECTED OUTPUT
#   A short paragraph in plain text comments at the end of the script.
#
# HINTS
#   • Lower SLA → conservative (slow-return) strategy: dense leaves, long
#     lifespan.
#   • Stress-tolerators and geophytes are typical low-SLA + threatened
#     species.
#   • Carmona et al. (2021) Nature show threatened plants tend to be slow-
#     pace.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group1_day1_*.png
# =============================================================================
