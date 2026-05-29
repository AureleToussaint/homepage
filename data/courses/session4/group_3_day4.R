# =============================================================================
# DAY 4  ·  GROUP 3  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q3 exercise for GROUP 3.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 5.
#
# RULES
#   • Use the ANNOTATED script (`day4_fish_landmarks_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group3_day4_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Compute FRic, FEve, FDis on 12 Amazonian fish communities, relate them to a degradation index, and test whether CWM(Body_elong) tracks current velocity.
#
# DATA / OBJECTS YOU NEED
#   • data/fish_traits_sp_means.csv
#   • data/fish_communities_amazon.csv
#   • data/fish_site_env.csv
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day4_fish_landmarks_annotated.R` (sections relevant to Q3) for context.
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
# Q3 — Functional diversity along the Amazonian degradation gradient
# =============================================================================
# -----------------------------------------------------------------------------
# Q3a — Compute FD indices for the 12 communities
# -----------------------------------------------------------------------------
# TASK
#   Run `dbFD()` on the scaled trait matrix and the community matrix. Identify
#   the site with highest and lowest functional diversity for each index
#   (FRic, FEve, FDiv, FDis).
#
# EXPECTED OUTPUT
#   Tibble of 12 rows × 6 columns (site, SR, FRic, FEve, FDiv, FDis).
#   Top/bottom sites identified.
#
# HINTS
#   • Make sure to align species columns of comm_mat with rows of trait_mat
#     before calling `dbFD()`.
#   • Set `m = 4` (PCoA axes) and `stand.x = FALSE` (already scaled).

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q3b — FD vs degradation
# -----------------------------------------------------------------------------
# TASK
#   Plot FRic and FDis as a function of degradation index
#   (`fish_site_env$degradation`). Fit linear models. Which index responds
#   more strongly? Why might they respond differently?
#
# EXPECTED OUTPUT
#   Two-panel scatter outputs/group3_day4_Q3b_FD_degradation.png with linear
#   fits and adjusted R². Discussion: FRic typically drops more abruptly (loss
#   of extreme species).
#
# HINTS
#   • `patchwork` to combine the two panels.
#   • Report `summary(lm)` adjusted R² and coefficient p-value.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q3c — CWM(Body_elong) vs current velocity
# -----------------------------------------------------------------------------
# TASK
#   Compute CWM(Body_elong) per site. Plot vs `current_velocity`. Test the
#   relationship with Spearman correlation (degradation and current may be
#   confounded).
#
# EXPECTED OUTPUT
#   PNG + Spearman r and p. Positive correlation expected: faster currents
#   favour elongated bodies (Villéger et al. 2010 Ecol. Appl.).
#
# HINTS
#   • Use `functcomp(trait_mat_sc, comm_mat)` and pull the Body_elong column.
#   • Spearman is robust to non-linearity.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group3_day4_*.png
# =============================================================================
