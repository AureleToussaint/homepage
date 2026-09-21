# =============================================================================
# DAY 4  ·  GROUP 1  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q1 exercise for GROUP 1.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 5.
#
# RULES
#   • Use the ANNOTATED script (`day4_fish_landmarks_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group1_day4_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Quantify measurement reliability for landmark-based traits (TEM, CV, ICC), identify which landmark is most error-prone and propose a biological reason.
#
# DATA / OBJECTS YOU NEED
#   • data/fish_landmarks_day4.csv
#   • irr package (`install.packages('irr')`)
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day4_fish_landmarks_annotated.R` (sections relevant to Q1) for context.
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
# Q1 — Measurement quality control: TEM, ICC and biological interpretation
# =============================================================================
# -----------------------------------------------------------------------------
# Q1a — Re-measure: TEM and CV per measurement
# -----------------------------------------------------------------------------
# TASK
#   Take 3 individuals from `fish_landmarks_day4.csv`. SIMULATE a second
#   measurement by adding ~5 % noise (in real practice this would be an actual
#   remeasurement). Compute the Technical Error of Measurement (TEM) and the
#   Coefficient of Variation (CV) per raw measurement (TL, BD, ED, ...).
#
# EXPECTED OUTPUT
#   Tibble with one row per measurement and columns mean, sd, TEM, CV.
#   Highlight the measurement with the highest CV.
#
# HINTS
#   • TEM = sqrt( sum( (M1 - M2)^2 ) / (2 * n) ).
#   • CV = sd / mean × 100 (in %).
#   • Use `mutate(across(TL_mm:BD_CP_mm, list(noisy = ~ . * (1 + rnorm(n(), 0,
#     0.05)))))` for the simulated remeasurement.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q1b — Intra-Class Correlation Coefficient (ICC)
# -----------------------------------------------------------------------------
# TASK
#   Use `irr::icc()` on the matrix of (M1, M2) per measurement. Compute the
#   ICC (two-way, agreement). Is the protocol reliable (ICC > 0.90)?
#
# EXPECTED OUTPUT
#   ICC value, 95 % CI, p-value. Comment on which measurements fall below the
#   0.90 threshold.
#
# HINTS
#   • `icc(cbind(M1, M2), model = 'twoway', type = 'agreement', unit =
#     'single')`.
#   • ICC > 0.90 = excellent; 0.75-0.90 = good; < 0.75 = unreliable.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q1c — Biological explanation of the highest-error measurement
# -----------------------------------------------------------------------------
# TASK
#   For the measurement with the highest CV from Q1a (typically MH or PA),
#   write a biological explanation: is it soft tissue? subjective endpoint?
#   age-dependent?
#
# EXPECTED OUTPUT
#   Comment block (4-6 lines) with a published reference.
#
# HINTS
#   • Mouth height (MH) = soft tissue, deformable; landmark 7 (dorsal fin
#     origin) is variable across taxa.
#   • Pectoral fin area (PA) requires outline tracing — much more subjective
#     than two-point distances.
#   • Cite Klingenberg (2016) on landmark digitisation error.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group1_day4_*.png
# =============================================================================
