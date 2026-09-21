# =============================================================================
# DAY 3 — Trait–Environment Relationships & Extinction Consequences
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil
# Dataset: 30 Amazonian fish species × 12 sites along a degradation gradient
# =============================================================================
# OBJECTIVES
#   Part 1 — Trait–Environment Relationships
#     - Build the L-Q-R three-table framework for a Neotropical fish dataset
#     - Run RLQ ordination and interpret biplots
#     - Perform fourth-corner permutation tests (models 2 & 4 + max-test)
#     - Visualise trait × environment relationships as FDR-corrected heatmaps
#
#   Part 2 — Extinction Consequences
#     - Simulate three extinction scenarios on Amazonian fish communities
#     - Track FRic and FDis loss along extinction sequences
#     - Compute FUSE scores (Functionally Unique, Specialised and Endangered)
#
# DATA FILES USED
#   data/fish_traits_sp_means.csv      — 30 species × 10 morphological traits
#   data/fish_communities_amazon.csv   — 12 sites × 30 species (abundances)
#   data/fish_site_env.csv             — 12 sites × 8 env. variables
#   data/fish_iucn_status.csv          — IUCN category per species
#
# PACKAGES REQUIRED
#   tidyverse, ade4, adespatial, FD, funrar, ggplot2, patchwork, ggrepel, corrplot
#
# KEY REFERENCES
#   Dray & Legendre (2008) Ecology 89: 3400–3412
#   Ter Braak et al. (2012) Methods Ecol. Evol. 3: 217–226
#   Mouillot et al. (2013) PLoS Biology 11: e1001569
#   Pimiento et al. (2020) Science Advances 6: eaay7650 — FUSE
#   Toussaint et al. (2016) Sci. Reports 6: 22125 — Neotropical fish FD
# =============================================================================

library(tidyverse)
library(ade4)        # dudi.coa(), dudi.pca(), rlq(), fourthcorner()
library(FD)          # dbFD()
library(funrar)      # uniqueness(), distinctiveness()
library(ggrepel)
library(patchwork)

theme_set(
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#E1F5EE"),
          strip.text       = element_text(colour = "#085041", face = "bold"))
)

if (!dir.exists("outputs")) dir.create("outputs")

# =============================================================================
# 1. LOAD AND PREPARE THE THREE TABLES
# =============================================================================
# In the RLQ framework:
#   R = sites × environmental variables   (12 × 8)
#   L = sites × species abundances        (12 × 30)
#   Q = species × functional traits       (30 × 10)
# Site names must match between R and L; species names between L and Q.

fish_traits <- read_csv("data/fish_traits_sp_means.csv", show_col_types = FALSE)
fish_comm   <- read_csv("data/fish_communities_amazon.csv", show_col_types = FALSE)
site_env    <- read_csv("data/fish_site_env.csv", show_col_types = FALSE)
iucn        <- read_delim("data/fish_iucn_status.csv",delim=";", show_col_types = FALSE)

glimpse(fish_traits)
glimpse(site_env)

# 1a. Table R — environmental variables (continuous, ecologically meaningful)
R_tab <- site_env |>
  select(site_id, current_velocity, water_depth, canopy_cover,
         pH, turbidity, conductivity, DO_mg_L, degradation) |>
  column_to_rownames("site_id")

# 1b. Table L — sites × species abundance matrix
L_tab <- fish_comm |>
  column_to_rownames("site_id")

# Sanity check: site rows must match
stopifnot(identical(rownames(R_tab), rownames(L_tab)))

# 1c. Table Q — species × functional traits (the 10 Villéger 2017 traits)
Q_tab <- fish_traits |>
  column_to_rownames("species")

# Species columns of L must match rows of Q (and in the SAME order for ade4)
stopifnot(setequal(colnames(L_tab), rownames(Q_tab)))
Q_tab <- Q_tab[colnames(L_tab), , drop = FALSE]

# =============================================================================
# 2. PRELIMINARY ORDINATIONS (REQUIRED BY RLQ)
# =============================================================================
# RLQ chains three weighted ordinations sharing the same row/column weights:
#   - Correspondence Analysis on L (sites × species)
#   - Weighted PCA on R using site weights from CA(L)
#   - Weighted PCA on Q using species weights from CA(L)

dudiL <- dudi.coa(L_tab, scannf = FALSE, nf = 4)

dudiR <- dudi.pca(R_tab, row.w = dudiL$lw,
                  scannf = FALSE, nf = 4, scale = TRUE)

dudiQ <- dudi.pca(Q_tab, row.w = dudiL$cw,
                  scannf = FALSE, nf = 4, scale = TRUE)

# =============================================================================
# 3. RLQ ANALYSIS
# =============================================================================
rlq_res <- rlq(dudiR, dudiL, dudiQ, scannf = FALSE, nf = 2)
summary(rlq_res)

# Variance explained by each RLQ axis
rlq_var <- rlq_res$eig / sum(rlq_res$eig)
cat("Variance explained per RLQ axis (%):\n")
print(round(rlq_var * 100, 1))

# Save for the exercises
saveRDS(rlq_res, "outputs/day3_rlq_result.rds")

# =============================================================================
# 4. RLQ BIPLOTS — sites, species, env. and trait loadings
# =============================================================================
# RLQ outputs (ade4 names):
#   $lR  : site scores       $lQ  : species scores
#   $l1  : env. loadings     $c1  : trait loadings

site_sc <- as.data.frame(rlq_res$lR) |>
  rownames_to_column("site_id") |>
  rename(Axis1 = 2, Axis2 = 3) |>
  left_join(site_env, by = "site_id")

sp_sc <- as.data.frame(rlq_res$lQ) |>
  rownames_to_column("species") |>
  rename(Axis1 = 2, Axis2 = 3)

env_load <- as.data.frame(rlq_res$l1) |>
  rownames_to_column("env_var") |>
  rename(load1 = 2, load2 = 3)

trait_load <- as.data.frame(rlq_res$c1) |>
  rownames_to_column("trait") |>
  rename(load1 = 2, load2 = 3)

# Plot A — Sites + environmental loading vectors
gg_sites <- ggplot(site_sc, aes(Axis1, Axis2)) +
  geom_hline(yintercept = 0, lty = 2, colour = "grey60") +
  geom_vline(xintercept = 0, lty = 2, colour = "grey60") +
  geom_point(aes(colour = degradation), size = 4) +
  scale_colour_viridis_c(option = "plasma", name = "Degradation") +
  geom_segment(data = env_load,
               aes(x = 0, y = 0, xend = load1 * 3, yend = load2 * 3),
               arrow = arrow(length = unit(0.2, "cm")),
               colour = "#085041", inherit.aes = FALSE) +
  geom_text_repel(data = env_load,
                  aes(x = load1 * 3, y = load2 * 3, label = env_var),
                  inherit.aes = FALSE, size = 3.2, colour = "#085041") +
  geom_text_repel(aes(label = site_id), size = 3) +
  labs(title = "RLQ — Sites in trait–environment space",
       x = paste0("Axis 1 (", round(rlq_var[1]*100, 1), "%)"),
       y = paste0("Axis 2 (", round(rlq_var[2]*100, 1), "%)"))

# Plot B — Species + trait loadings
gg_sp <- ggplot(sp_sc, aes(Axis1, Axis2)) +
  geom_hline(yintercept = 0, lty = 2, colour = "grey60") +
  geom_vline(xintercept = 0, lty = 2, colour = "grey60") +
  geom_point(size = 2.5, colour = "#0A4A5A", alpha = 0.85) +
  geom_segment(data = trait_load,
               aes(x = 0, y = 0, xend = load1 * 3, yend = load2 * 3),
               arrow = arrow(length = unit(0.2, "cm")),
               colour = "#B23A48", inherit.aes = FALSE) +
  geom_text_repel(data = trait_load,
                  aes(x = load1 * 3, y = load2 * 3, label = trait),
                  size = 3.2, colour = "#B23A48", inherit.aes = FALSE) +
  geom_text_repel(aes(label = species), size = 2.6, max.overlaps = 20) +
  labs(title = "RLQ — Species in trait–environment space",
       x = paste0("Axis 1 (", round(rlq_var[1]*100, 1), "%)"),
       y = paste0("Axis 2 (", round(rlq_var[2]*100, 1), "%)"))

p_rlq <- gg_sites + gg_sp + plot_layout(ncol = 2)
ggsave("outputs/day3_rlq_biplots.png", p_rlq,
       width = 13, height = 6, dpi = 300)

# =============================================================================
# 5. FOURTH-CORNER TEST — significance of trait × environment associations
# =============================================================================
# Why two permutation models?
#   - model 2: permute rows of L     -> tests trait response to environment
#   - model 4: permute columns of L  -> tests environment response to traits
# Used alone, each one INFLATES Type-I error. The Ter Braak et al. (2012)
# "max-test" combines them (takes the MAX of the two p-values cell-by-cell)
# -> controlled Type-I error rate.
#
# Recommended workflow in ade4: run BOTH models, then combine with
# combine.4thcorner() — this returns a 4thcorner object whose $tabD2 is a
# clean data.frame (one row per env x trait pair), ready for ggplot.
#
# WARNING about other functions:
#   * fourthcorner2(R, L, Q): a DIFFERENT function — computes a single global
#     Sxy statistic. It takes data frames, not two fourthcorner results.
#   * fourthcorner(..., modeltype = 6): performs the max-test in one call
#     but returns a krandtest in $tabD2 (NOT a data.frame), which complicates
#     downstream tidyverse workflows. We avoid it here.

set.seed(42)
fc_m2 <- fourthcorner(R_tab, L_tab, Q_tab,
                      modeltype = 2, nrepet = 999,
                      p.adjust.method.G = "fdr",
                      p.adjust.method.D = "fdr")

set.seed(42)
fc_m4 <- fourthcorner(R_tab, L_tab, Q_tab,
                      modeltype = 4, nrepet = 999,
                      p.adjust.method.G = "fdr",
                      p.adjust.method.D = "fdr")

# Combined max-test (returns a 4thcorner with $tabD2 as data.frame)
fc_combined <- combine.4thcorner(fc_m2, fc_m4)

saveRDS(fc_m2,       "outputs/day3_fc_m2.rds")
saveRDS(fc_m4,       "outputs/day3_fc_m4.rds")
saveRDS(fc_combined, "outputs/day3_fourthcorner_combined.rds")

# =============================================================================
# 6. HEATMAP OF TRAIT × ENVIRONMENT ASSOCIATIONS
# =============================================================================
# fc_combined$tabD2 contains one row per (env, trait) cell.
# In recent ade4, several of its columns can be LIST-COLUMNS where each cell
# is a 1-element list, sometimes empty. We use a single bulletproof helper
# (`pull_vec`) that always returns a clean atomic vector of the expected
# length, regardless of whether the column is atomic, list-of-1, list-of-N,
# or contains NULLs.

n_traits <- ncol(Q_tab)
n_envs   <- ncol(R_tab)
expected_n <- n_envs * n_traits

# Bulletproof column flattener
pull_vec <- function(x, expected_n, mode = c("numeric", "character")) {
  mode <- match.arg(mode)
  na_value <- if (mode == "numeric") NA_real_ else NA_character_
  caster   <- if (mode == "numeric") as.numeric else as.character

  # Try the obvious path first
  out <- tryCatch(
    suppressWarnings(caster(unlist(x, use.names = FALSE))),
    error = function(e) NULL
  )
  if (!is.null(out) && length(out) == expected_n) return(out)

  # Fallback: iterate element-wise, returning NA for empty/NULL entries
  out <- vapply(seq_len(expected_n), function(i) {
    e <- if (i <= length(x)) x[[i]] else NULL
    if (is.null(e) || length(e) == 0) na_value
    else suppressWarnings(caster(unlist(e, use.names = FALSE))[1])
  }, FUN.VALUE = if (mode == "numeric") numeric(1) else character(1))
  out
}

tabD2 <- fc_combined$tabD2
obs_vec  <- pull_vec(tabD2$obs,        expected_n, "numeric")
padj_vec <- pull_vec(tabD2$adj.pvalue, expected_n, "numeric")

# We don't parse $names at all — we rely on the ade4 convention that tabD2
# rows are ordered env-major / trait-minor:
#   row 1: env_1 × trait_1
#   row 2: env_1 × trait_2
#   ...
#   row n_traits + 1: env_2 × trait_1, etc.
fc_tab <- tibble(
  env   = rep(colnames(R_tab), each  = n_traits),
  trait = rep(colnames(Q_tab), times = n_envs),
  stat  = obs_vec,
  p.adj = padj_vec
) |>
  mutate(sign = case_when(p.adj < 0.001 ~ "***",
                          p.adj < 0.01  ~ "**",
                          p.adj < 0.05  ~ "*",
                          TRUE          ~ ""),
         fill_r = ifelse(p.adj < 0.05, stat, NA_real_))

# Quick sanity check (warning only, does not abort the script)
if (all(is.na(fc_tab$stat))) {
  warning("All test statistics are NA — check fc_combined$tabD2$obs structure ",
          "with str(fc_combined$tabD2, max.level = 2).")
}

p_heatmap <- ggplot(fc_tab, aes(env, trait, fill = fill_r)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sign), size = 4) +
  scale_fill_gradient2(low  = "#2C7BB6", mid = "white", high = "#D7191C",
                       midpoint = 0, na.value = "grey92",
                       name = "r (FDR-sig.)") +
  labs(title    = "Fourth-corner — trait × environment (FDR-corrected, max-test)",
       subtitle = "* p<0.05  ** p<0.01  *** p<0.001",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("outputs/day3_fourthcorner_heatmap.png", p_heatmap,
       width = 8, height = 5.5, dpi = 300)

saveRDS(fc_tab, "outputs/day3_fc_heatmap_data.rds")

# =============================================================================
# 7. EXTINCTION SCENARIOS — functional consequences of species loss
# =============================================================================
# Compute baseline functional diversity, then iteratively remove species
# according to three biologically motivated scenarios (Mouillot et al. 2013).

# Trait matrix scaled (different units across the 10 traits)
trait_mat <- fish_traits |>
  column_to_rownames("species") |>
  as.matrix()
trait_mat_sc <- scale(trait_mat)

# Community matrix (sites × species)
comm_mat <- as.matrix(L_tab)

# Baseline FD on the full community
fd_full <- dbFD(x = trait_mat_sc, a = comm_mat,
                calc.FRic = TRUE, calc.FDiv = TRUE, calc.FGR = FALSE,
                stand.x = FALSE, m = 4, print.pco = FALSE,
                messages = FALSE)

baseline <- tibble(site_id   = names(fd_full$nbsp),
                   FRic_full = fd_full$FRic,
                   FDis_full = fd_full$FDis)

# Build the three extinction sequences
abund_total <- colSums(comm_mat)

iucn_w <- iucn |>
  mutate(weight = recode(iucn_status,
                         "LC" = 0.0,  "NT" = 0.10,
                         "VU" = 0.40, "EN" = 0.667,
                         "CR" = 0.999, "DD" = 0.30)) |>
  arrange(desc(weight))

set.seed(20260511)
seq_random  <- sample(rownames(trait_mat_sc))
seq_lowabd  <- names(sort(abund_total, decreasing = FALSE))
seq_threat  <- iucn_w$species

# Function: compute FD along an extinction sequence
# Defensive: catches dbFD() failures (e.g. when the remaining trait matrix
# becomes too small for the requested number of PCoA axes m).
fd_along <- function(seq_order, trait_m, comm_m){
  res <- vector("list", length(seq_order))
  for (i in seq_along(seq_order)) {
    keep <- setdiff(rownames(trait_m), seq_order[seq_len(i - 1)])
    if (length(keep) < 4) break
    cm <- comm_m[, keep, drop = FALSE]
    cm <- cm[rowSums(cm) > 0, , drop = FALSE]
    if (nrow(cm) < 1) break
    # Adapt m to the available number of dimensions
    m_use <- min(4, length(keep) - 1)
    fd <- tryCatch(
      dbFD(trait_m[keep, ], cm,
           calc.FRic = TRUE, stand.x = FALSE, m = m_use,
           messages = FALSE, print.pco = FALSE),
      error = function(e) NULL
    )
    if (is.null(fd)) next  # skip this iteration, continue removing species
    res[[i]] <- tibble(step = i - 1L,
                       species_lost = if (i == 1L) NA_character_ else seq_order[i - 1L],
                       FRic_mean = mean(fd$FRic, na.rm = TRUE),
                       FDis_mean = mean(fd$FDis, na.rm = TRUE),
                       SR_mean   = mean(fd$nbsp,  na.rm = TRUE))
  }
  bind_rows(res)
}

ext_random <- fd_along(seq_random, trait_mat_sc, comm_mat) |>
  mutate(scenario = "A · Random")
ext_lowabd <- fd_along(seq_lowabd, trait_mat_sc, comm_mat) |>
  mutate(scenario = "B · Least abundant first")
ext_threat <- fd_along(seq_threat, trait_mat_sc, comm_mat) |>
  mutate(scenario = "C · Most threatened first")

ext_all <- bind_rows(ext_random, ext_lowabd, ext_threat)
saveRDS(ext_all, "outputs/day3_extinction_curves.rds")

p_ext <- ggplot(ext_all, aes(step, FRic_mean, colour = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = c("#7FB069", "#F4A261", "#B23A48")) +
  labs(title    = "Functional space erosion under three extinction scenarios",
       subtitle = "Amazonian fish — 30 species, 12 sites",
       x = "Number of species removed",
       y = "Mean FRic across sites",
       colour = NULL) +
  theme(legend.position = "bottom")

ggsave("outputs/day3_extinction_curves.png", p_ext,
       width = 8, height = 5, dpi = 300)

# =============================================================================
# 8. FUSE — Functionally Unique, Specialised and Endangered
# =============================================================================
# FUSE = FUn_norm × FSpec_norm × IUCN_weight   (Pimiento et al. 2020)
#   FUn   = nearest-neighbour distance in trait space
#   FSpec = distance from the centroid of the global functional space
#   IUCN_weight (Mooers et al. 2008): LC=0, NT=0.1, VU=0.4, EN=0.667, CR=0.999

trait_dist <- dist(trait_mat_sc)

trait_dm <- as.matrix(trait_dist)
diag(trait_dm) <- NA
FUn <- apply(trait_dm, 1, min, na.rm = TRUE)

centroid <- colMeans(trait_mat_sc)
FSpec    <- apply(trait_mat_sc, 1, function(x) sqrt(sum((x - centroid)^2)))

norm01 <- function(x) (x - min(x, na.rm = TRUE)) /
                      (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

fuse_df <- tibble(species = rownames(trait_mat_sc),
                  FUn = FUn, FUn_n = norm01(FUn),
                  FSpec = FSpec, FSpec_n = norm01(FSpec)) |>
  left_join(iucn, by = "species") |>
  mutate(IUCN_w = recode(iucn_status,
                         "LC" = 0.0,  "NT" = 0.10,
                         "VU" = 0.40, "EN" = 0.667,
                         "CR" = 0.999, "DD" = 0.30),
         FUSE = FUn_n * FSpec_n * IUCN_w) |>
  arrange(desc(FUSE))

saveRDS(fuse_df,  "outputs/day3_fuse_scores.rds")
write_csv(fuse_df, "outputs/day3_fuse_scores.csv")

p_fuse <- ggplot(fuse_df, aes(FUn_n, FSpec_n)) +
  geom_point(aes(size = FUSE, colour = iucn_status), alpha = 0.85) +
  geom_text_repel(aes(label = species), size = 2.6, max.overlaps = 25) +
  scale_size(range = c(2, 9), name = "FUSE") +
  scale_colour_manual(values = c("LC" = "grey60", "NT" = "#F4A261",
                                 "VU" = "#E76F51", "EN" = "#9D2226",
                                 "CR" = "#3A0CA3", "DD" = "#5D5D5D"),
                      name = "IUCN") +
  labs(title    = "FUSE — Amazonian fish: functional uniqueness × specialisation",
       subtitle = "Bubble size proportional to FUSE = FUn × FSpec × IUCN_weight",
       x = "Functional Uniqueness (normalised)",
       y = "Functional Specialisation (normalised)")

ggsave("outputs/day3_fuse_scores.png", p_fuse,
       width = 9, height = 6, dpi = 300)

# =============================================================================
# 9. SAVE OBJECTS NEEDED BY THE EXERCISES
# =============================================================================
saveRDS(list(R = R_tab, L = L_tab, Q = Q_tab,
             trait_mat_sc = trait_mat_sc,
             comm_mat     = comm_mat,
             iucn         = iucn,
             baseline     = baseline,
             fuse_df      = fuse_df),
        "outputs/day3_objects.rds")

cat("\n=== DAY 3 — outputs written to outputs/ ===\n")
cat("Top-5 FUSE species (Amazonian fish):\n")
print(fuse_df |> select(species, iucn_status, FUSE) |> slice_head(n = 5))

# =============================================================================
# EXERCISES — see exercises_day3.R
#   Q1 RLQ interpretation     → re-run without Body_size, top/bottom 5 species, Moran's I
#   Q2 Fourth-corner critical → 3 strongest associations, Type-I inflation, linearity
#   Q3 Extinction scenarios   → 4th 'most common first', %FRic loss for CR+EN
#   Q4 FUSE in context        → top-5 priorities, regional weights, mechanistic discussion
# =============================================================================
