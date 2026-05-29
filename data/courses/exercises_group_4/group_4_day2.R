# =============================================================================
# DAY 2  ·  GROUP 4  —  Student Exercise Sheet
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil (May 2026)
# =============================================================================
# This script is the Q4 exercise for GROUP 4.
# Submit your annotated R script with figures saved to outputs/ before the morning of Day 3.
#
# RULES
#   • Use the ANNOTATED script (`day2_functional_diversity_annotated.R`) and dataset references below.
#   • Comment your code: explain WHY you choose each step, not only WHAT it does.
#   • Save every figure to outputs/group4_day2_<short>.png.
#   • Do NOT copy code from another group.
#
# LEARNING GOAL
#   Distinguish three related but distinct measures of functional uniqueness: global originality, nearest-neighbour FUn, and TPD-based rarity.
#
# DATA / OBJECTS YOU NEED
#   • trait_mat_sc, sp_traits, fs (funspace object)
#
# REFERENCE TO ANNOTATED SCRIPT
#   See `day2_functional_diversity_annotated.R` (sections relevant to Q4) for context.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. SETUP — load packages and rebuild required objects
# -----------------------------------------------------------------------------
library(tidyverse)
library(FD)
library(funspace)
library(vegan)
library(ggrepel)
library(ggforce)
library(patchwork)
library(corrplot)

theme_set(theme_bw(base_size = 12))
if (!dir.exists("outputs")) dir.create("outputs")

# Pre-requisite: outputs/day1_sp_traits_clean.rds must exist (Day 1 output).
sp_traits <- readRDS("outputs/day1_sp_traits_clean.rds")
env       <- read_csv("data/02_environment.csv",           show_col_types = FALSE)
comm_raw  <- read_csv("data/03_community_composition.csv", show_col_types = FALSE)

# Build matrices (same as in the annotated Day 2 script)
trait_mat <- sp_traits |>
  select(species, SLA, LDMC, LeafN, Height, SeedMass) |>
  mutate(across(c(SLA, LDMC, Height, SeedMass), log)) |>
  column_to_rownames("species") |> as.matrix()
trait_mat_sc <- scale(trait_mat)
comm_mat <- comm_raw |> column_to_rownames("plot_id") |> as.matrix()

# Common FD output (used by Q1-Q4)
fd_full <- dbFD(trait_mat_sc, comm_mat, calc.FRic = TRUE, calc.FDiv = TRUE,
                stand.x = FALSE, m = 3, print.pco = FALSE, messages = FALSE)

# =============================================================================
# Q4 — Functional rarity, redundancy and originality
# =============================================================================
# -----------------------------------------------------------------------------
# Q4a — Compute functional originality
# -----------------------------------------------------------------------------
# TASK
#   For each species, compute the MEAN distance to ALL other species in the
#   PCoA / scaled-trait space. Plot originality versus IUCN threat weight. Are
#   functionally original species disproportionately threatened?
#
# EXPECTED OUTPUT
#   Tibble of originality per species + scatter plot saved as
#   outputs/group4_day2_Q4a_originality.png. Positive correlation expected:
#   original species more likely threatened.
#
# HINTS
#   • `d <- as.matrix(dist(trait_mat_sc))`; originality = `rowMeans(d)` after
#     setting `diag(d) <- NA`.
#   • `cor.test(orig, iucn_w, method = 'spearman')`.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q4b — Map originality onto the funspace
# -----------------------------------------------------------------------------
# TASK
#   Use `funspaceGAM(fs, var = originality, family = gaussian())` to show
#   which regions of trait space host the most original species. Identify the
#   2 zones of highest originality.
#
# EXPECTED OUTPUT
#   GAM map outputs/group4_day2_Q4b_originality_map.png. Highest originality
#   typically at the periphery of the global cloud.
#
# HINTS
#   • Originality is high in regions far from the centroid and from any
#     neighbour — usually the corners of the LES space.

# YOUR CODE HERE
# ──────────────


# -----------------------------------------------------------------------------
# Q4c — Discussion: rarity vs originality vs FUn
# -----------------------------------------------------------------------------
# TASK
#   Write 100+ words distinguishing the three concepts. Use Day 3's FUn
#   (nearest-neighbour distance) as a contrast point. When would you prefer
#   one over the others for conservation prioritisation?
#
# EXPECTED OUTPUT
#   Comment block at the end of the script.
#
# HINTS
#   • Originality = mean distance to ALL others (global isolation).
#   • FUn = distance to the NEAREST neighbour (local isolation).
#   • Rarity (TPD-based) = density-weighted isolation.
#   • FUn is sensitive to outliers; originality dilutes them; TPD-rarity is
#     robust but requires a probabilistic framework.

# YOUR CODE HERE
# ──────────────


# =============================================================================
# END OF EXERCISE — submit your annotated script + outputs/group4_day2_*.png
# =============================================================================
