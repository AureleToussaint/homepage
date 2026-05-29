# =============================================================================
# DAY 2 — Community-Weighted Means and Functional Diversity Indices
# Course: Trait-Based Ecology — Theory and Applications in R
# =============================================================================
# OBJECTIVE
#   - Compute CWM from community composition and trait matrices
#   - Calculate FRic, FEve, FDiv, FDis using {FD}
#   - Visualise functional space and interpret FD indices
#   - Assess sensitivity of FD indices to trait selection and sample size
#
# DATA FILES USED
#   outputs/day1_sp_traits_clean.rds      (from Day 1)
#   data/03_community_composition.csv
#   data/02_environment.csv
#
# PACKAGES REQUIRED
#   tidyverse, FD, vegan, ggplot2, patchwork, ggforce
# =============================================================================

library(tidyverse)
library(FD)        # dbFD(), functcomp()
library(vegan)     # diversity(), decostand()
library(patchwork)
library(ggforce)   # geom_mark_ellipse()

theme_set(
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#E1F5EE"),
          strip.text       = element_text(colour = "#085041", face = "bold"))
)

# ── 1. Load data ───────────────────────────────────────────────────────────────

sp_traits <- readRDS("outputs/day1_sp_traits_clean.rds")
comm_raw  <- read_csv("data/03_community_composition.csv", show_col_types = FALSE)
env       <- read_csv("data/02_environment.csv",           show_col_types = FALSE)

# ── 2. Prepare matrices for {FD} ──────────────────────────────────────────────
# dbFD() expects:
#   x   = trait matrix   (species × traits, row names = species names)
#   a   = abundance matrix (sites × species, row names = site names)
# IMPORTANT: species names must match EXACTLY between x and a

# 2a. Trait matrix (log-transform skewed traits for Euclidean distance quality)
trait_mat <- sp_traits |>
  select(species, SLA, LDMC, LeafN, Height, SeedMass) |>
  mutate(
    SLA      = log(SLA),
    LDMC     = log(LDMC),
    Height   = log(Height),
    SeedMass = log(SeedMass)
  ) |>
  column_to_rownames("species") |>
  as.matrix()

# Scale traits to unit variance — mandatory when traits have different units
# Without scaling, Height (cm, large values) would dominate the functional space
trait_mat_sc <- scale(trait_mat)
cat("Trait matrix dimensions:", dim(trait_mat_sc), "\n")

# 2b. Community matrix (sites × species)
comm_mat <- comm_raw |>
  column_to_rownames("plot_id") |>
  as.matrix()

# Ensure species order is identical between trait and community matrices
shared_sp <- intersect(rownames(trait_mat_sc), colnames(comm_mat))
cat("Species in trait matrix:", nrow(trait_mat_sc), "\n")
cat("Species in community matrix:", ncol(comm_mat), "\n")
cat("Shared species:", length(shared_sp), "\n")

# Subset to shared species
trait_mat_sc <- trait_mat_sc[shared_sp, ]
comm_mat     <- comm_mat[, shared_sp]

# Remove plots with zero total cover (can happen after species filtering)
comm_mat <- comm_mat[rowSums(comm_mat) > 0, ]


# ── 3. Community-Weighted Mean (CWM) ──────────────────────────────────────────
# CWM_trait = Σ (relative abundance_i × trait_i) for all species i in the community
# Equivalent to the abundance-weighted average trait value
# Interpretation: reflects the dominant strategy in the community (mass-ratio hypothesis)

# Manual computation (pedagogical)
# Note: comm_mat rows already sum to 1 (relative cover), so no need to divide
cwm_manual <- comm_mat %*% trait_mat_sc
cwm_manual <- as.data.frame(cwm_manual)
cwm_manual$plot_id <- rownames(cwm_manual)

# Using functcomp() from {FD} — same result, more convenient
cwm_fd <- functcomp(
  x = trait_mat_sc,
  a = comm_mat,
  CWM.type = "all"   # returns mean for continuous traits
)

# Verify manual == package results
cat("Max difference between manual and FD CWM (SLA):",
    max(abs(cwm_manual$SLA - cwm_fd$SLA)), "\n")   # should be ≈ 0

# Join CWM with environmental data
cwm <- cwm_fd |>
  rownames_to_column("plot_id") |>
  left_join(env, by = "plot_id")


# ── 4. CWM along the elevational gradient ────────────────────────────────────
# Prediction: CWM_SLA should decrease with elevation (conservative strategy prevails)
# Prediction: CWM_LeafN should decrease with elevation (nutrient limitation)

p_cwm_sla <- ggplot(cwm, aes(x = elevation, y = SLA)) +
  geom_point(aes(colour = disturbance), size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, colour = "#085041") +
  scale_colour_gradient(low = "#E1F5EE", high = "#D85A30",
                        name = "Disturbance") +
  labs(x = "Elevation (m)", y = "CWM log(SLA) — scaled",
       title = "A — CWM(SLA) along elevational gradient") +
  annotate("text", x = 2600, y = max(cwm$SLA)*0.9,
           label = paste0("r = ", round(cor(cwm$elevation, cwm$SLA), 2)),
           size = 4, colour = "#085041")

p_cwm_leafn <- ggplot(cwm, aes(x = elevation, y = LeafN)) +
  geom_point(aes(colour = disturbance), size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, colour = "#085041") +
  scale_colour_gradient(low = "#E1F5EE", high = "#D85A30",
                        name = "Disturbance") +
  labs(x = "Elevation (m)", y = "CWM log(LeafN) — scaled",
       title = "B — CWM(Leaf N) along elevational gradient")

p_cwm_sla + p_cwm_leafn +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


# ── 5. Computing FD indices with dbFD() ───────────────────────────────────────
# dbFD() computes all standard FD indices from a trait matrix and abundance matrix
# Key arguments:
#   calc.FRic  — compute functional richness (convex hull volume)? TRUE/FALSE
#   calc.FDiv  — compute FDiv and FEve? TRUE/FALSE
#   m          — number of PCoA axes to retain (default: all non-negative)
#   stand.x    — standardise trait matrix? (we already did it manually)
#   w          — trait weights (default: equal weights)

cat("Computing FD indices — this may take 1–2 minutes...\n")
fd_results <- dbFD(
  x         = trait_mat_sc,
  a         = comm_mat,
  calc.FRic = TRUE,
  calc.FDiv = TRUE,
  stand.x   = FALSE,   # already scaled
  m         = 3        # retain 3 PCoA axes (sufficient for 5 traits)
)

# fd_results is a list — extract key indices
fd_df <- data.frame(
  plot_id = rownames(comm_mat),
  SR      = fd_results$nbsp,      # species richness
  FRic    = fd_results$FRic,      # functional richness (convex hull volume)
  FEve    = fd_results$FEve,      # functional evenness
  FDiv    = fd_results$FDiv,      # functional divergence
  FDis    = fd_results$FDis,      # functional dispersion (=weighted mean dist to centroid)
  RaoQ    = fd_results$RaoQ       # Rao's quadratic entropy
) |>
  left_join(env, by = "plot_id")

cat("FD indices computed for", nrow(fd_df), "plots.\n")

# Quick summary statistics
fd_df |>
  select(SR, FRic, FEve, FDiv, FDis) |>
  summary()


# ── 6. Relationships between FD indices and elevation ─────────────────────────
# Expected: FRic peaks at intermediate elevations (hump-shaped pattern)?
# FDis decreases with elevation (filtering increases trait convergence)?

p_fric <- ggplot(fd_df, aes(x = elevation, y = FRic)) +
  geom_point(aes(colour = SR), size = 2.5) +
  geom_smooth(method = "loess", se = TRUE, colour = "#185FA5") +
  scale_colour_viridis_c(name = "Species richness") +
  labs(x = "Elevation (m)", y = "FRic",
       title = "A — Functional richness along elevation")

p_fdis <- ggplot(fd_df, aes(x = elevation, y = FDis)) +
  geom_point(aes(colour = disturbance), size = 2.5) +
  geom_smooth(method = "loess", se = TRUE, colour = "#185FA5") +
  scale_colour_gradient(low = "#E1F5EE", high = "#D85A30",
                        name = "Disturbance") +
  labs(x = "Elevation (m)", y = "FDis",
       title = "B — Functional dispersion along elevation")

p_feve <- ggplot(fd_df, aes(x = FEve, y = FDiv)) +
  geom_point(aes(colour = elevation), size = 2.5) +
  scale_colour_viridis_c(name = "Elevation (m)", option = "C") +
  labs(x = "Functional evenness (FEve)", y = "Functional divergence (FDiv)",
       title = "C — FEve vs. FDiv relationship") +
  annotate("text", x = min(fd_df$FEve, na.rm=TRUE) + 0.05,
           y = max(fd_df$FDiv, na.rm=TRUE),
           label = paste0("r = ", round(cor(fd_df$FEve, fd_df$FDiv,
                                            use="complete.obs"), 2)),
           size = 4)

(p_fric + p_fdis) / p_feve


# ── 7. Visualising functional space ───────────────────────────────────────────
# We use the first 2 PCoA axes from the trait distance matrix
# to visualise the functional space occupied by communities

# PCoA of trait matrix (Euclidean distance after scaling)
trait_dist <- dist(trait_mat_sc, method = "euclidean")
pcoa_res   <- cmdscale(trait_dist, k = 2, eig = TRUE)

# % variance explained by each axis
eig_pct <- round(100 * pcoa_res$eig[1:2] / sum(pcoa_res$eig[pcoa_res$eig > 0]), 1)
cat(sprintf("PCoA Axis 1: %.1f%% | Axis 2: %.1f%%\n", eig_pct[1], eig_pct[2]))

# Species scores in functional space
sp_pcoa <- as.data.frame(pcoa_res$points)
colnames(sp_pcoa) <- c("PC1", "PC2")
sp_pcoa <- sp_pcoa |>
  rownames_to_column("species") |>
  left_join(sp_traits |> select(species, PFT_label), by = "species")

# Plot: species in functional space, coloured by PFT
ggplot(sp_pcoa, aes(x = PC1, y = PC2, colour = PFT_label, label = species)) +
  ggforce::geom_mark_ellipse(
    aes(fill = PFT_label, colour = PFT_label),
    alpha = 0.12, linewidth = 0.6
  ) +
  geom_point(size = 3.5) +
  ggrepel::geom_text_repel(size = 2.6, max.overlaps = 12,
                           segment.colour = "grey70") +
  scale_colour_brewer(palette = "Set2", name = "PFT",
                      aesthetics = c("colour", "fill")) +
  labs(
    x     = paste0("PCoA 1 (", eig_pct[1], "% variance)"),
    y     = paste0("PCoA 2 (", eig_pct[2], "% variance)"),
    title = "Functional space — 30 alpine grassland species",
    subtitle = "Ellipses = 95% confidence regions per PFT"
  )


# ── 8. Sensitivity analysis — effect of trait number ─────────────────────────
# Key concern: does FRic change substantially when we add/remove traits?
# This is a major source of methodological uncertainty (Mouillot et al. 2021)

trait_combos <- list(
  "2 traits (SLA + LDMC)"            = c("SLA", "LDMC"),
  "3 traits (+ LeafN)"               = c("SLA", "LDMC", "LeafN"),
  "4 traits (+ Height)"              = c("SLA", "LDMC", "LeafN", "Height"),
  "5 traits (+ SeedMass)"            = c("SLA", "LDMC", "LeafN", "Height", "SeedMass")
)

fdric_sensitivity <- map_dfr(names(trait_combos), function(combo_name) {
  selected <- trait_combos[[combo_name]]
  mat_sub  <- trait_mat_sc[, selected, drop = FALSE]

  fd_sub <- suppressMessages(
    dbFD(x = mat_sub, a = comm_mat, calc.FRic = TRUE,
         calc.FDiv = FALSE, stand.x = FALSE, m = min(2, ncol(mat_sub)))
  )
  data.frame(
    combo   = combo_name,
    plot_id = names(fd_sub$FRic),
    FRic    = fd_sub$FRic
  )
})

# Visualise variance in FRic across trait combinations
ggplot(fdric_sensitivity, aes(x = combo, y = FRic, fill = combo)) +
  geom_violin(alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = 21, colour = "grey30") +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  labs(
    x     = NULL, y = "Functional richness (FRic)",
    title = "Sensitivity of FRic to number of traits included",
    caption = "Each violin = distribution of FRic across 40 plots."
  ) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# CRITICAL INTERPRETATION: FRic values are not comparable across different
# numbers of traits. Always use the same trait set across all communities.


# ── 9. Species richness vs. functional diversity ──────────────────────────────
# A key question: is FD simply a proxy for SR, or does it carry independent info?

ggplot(fd_df, aes(x = SR, y = FDis)) +
  geom_point(aes(colour = elevation), size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "#185FA5", linewidth = 1) +
  scale_colour_viridis_c(name = "Elevation (m)") +
  labs(
    x     = "Species richness (SR)",
    y     = "Functional dispersion (FDis)",
    title = "SR vs. FDis — are they independent?",
    caption = paste0("r = ", round(cor(fd_df$SR, fd_df$FDis,
                                       use = "complete.obs"), 2))
  )
# If r is high (> 0.7), consider using the FD residuals from SR
# (i.e., functional diversity unexplained by species richness)


# ── 10. Save outputs ─────────────────────────────────────────────────────────
saveRDS(fd_df,   "outputs/day2_fd_indices.rds")
saveRDS(cwm,     "outputs/day2_cwm.rds")
saveRDS(sp_pcoa, "outputs/day2_pcoa_scores.rds")
cat("Saved: outputs/day2_*.rds\n")

# =============================================================================
# EXERCISE — Students complete the following
# =============================================================================
# Q1. Compute CWM for plant height. Plot CWM(Height) as a function of
#     disturbance index. Interpret: which PFT drives the pattern?
#
# Q2. Compute FRic using ONLY the two LES axes (SLA and LDMC).
#     Compare communities with high vs. low FEve. What does low FEve
#     indicate about species packing in functional space?
#
# Q3. FRic is sensitive to the number of species in the community.
#     Compute the Pearson correlation between SR and FRic across plots.
#     Propose a method to separate the effect of SR from FRic.
#     (Hint: look at ?dbFD argument 'calc.FRic' and Villeger et al. 2008)
#
# Q4. Identify the three plots with the highest FDis and the three with the
#     lowest FDis. What are their elevations and disturbance levels?
#     Formulate a biological hypothesis to explain this pattern.
# =============================================================================
