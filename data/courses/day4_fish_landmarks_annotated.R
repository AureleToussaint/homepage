# =============================================================================
# DAY 4 — Fish Morphological Traits & Functional Diversity
# Course: Trait-Based Approach in Community Ecology — Belém, Brazil
# Dataset: 12 Neotropical fish species (raw landmarks) + 12 sites × 30 species
# =============================================================================
# OBJECTIVES
#   - Derive the 10 standardised morphological traits (Villéger et al. 2017)
#     from raw landmark/anatomical measurements
#   - Aggregate individual measurements to species means and inspect ITV
#   - Compute functional diversity indices (FRic, FEve, FDiv, FDis) on
#     Amazonian fish communities
#   - Relate FD to environmental degradation along the gradient
#   - Visualise the Amazonian fish functional space and trophic guilds
#
# DATA FILES USED
#   data/fish_landmarks_day4.csv      — raw measurements per individual
#                                       (12 spp × 5 ind. = 60 fish)
#   data/fish_traits_sp_means.csv     — pre-computed species means (30 spp)
#   data/fish_communities_amazon.csv  — sites × species matrix
#   data/fish_site_env.csv            — site environmental metadata
#   data/fish_species_metadata.csv    — common name, family, trophic guild
#
# PACKAGES
#   tidyverse, FD, vegan, ggrepel, patchwork, corrplot, irr
#
# KEY REFERENCES
#   Villéger et al. (2017) Ecography 40: 1–7 — fish trait protocol
#   Brosse et al. (2021) Sci. Data 8: 254 — FishMorph
#   Toussaint et al. (2016) Sci. Reports 6: 22125 — Neotropical FD
# =============================================================================

library(tidyverse)
library(FD)            # dbFD()
library(vegan)         # rda(), procrustes()
library(ggrepel)
library(patchwork)
library(corrplot)

theme_set(
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#E0F0F5"),
          strip.text       = element_text(colour = "#0A4A5A", face = "bold"))
)

if (!dir.exists("outputs")) dir.create("outputs")

# =============================================================================
# 1. LOAD RAW LANDMARK MEASUREMENTS
# =============================================================================
# Output of the ImageJ practical: one row per individual, columns are
# anatomical distances in millimetres.

raw <- read_csv("data/fish_landmarks_day4.csv", show_col_types = FALSE)
glimpse(raw)
# 60 specimens × measurements (TL, BD, BW, ED, MH, MW, HMd, PH, PA, CFL, CFA,
#                              CPL, PFP, BD_CP)

# Quick QA — visualise the distribution of body lengths per species
p_qa <- ggplot(raw, aes(species, TL_mm)) +
  geom_boxplot(fill = "#A8D5BA", outlier.colour = "#B23A48") +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
  coord_flip() +
  labs(title = "Body size (TL) variation across measured Amazonian species",
       y = "Total length (mm)", x = NULL)
ggsave("outputs/day4_qa_body_size.png", p_qa, width = 8, height = 5, dpi = 300)

# =============================================================================
# 2. COMPUTE THE 10 MORPHOLOGICAL TRAITS (VILLÉGER ET AL. 2017)
# =============================================================================
# All traits are dimensionless ratios — comparable across species of very
# different body sizes (a 3 cm Astyanax to a 3 m Arapaima).

fish_traits_ind <- raw |>
  mutate(
    # ── Locomotion (6 traits) ─────────────────────────────────────────
    Body_size   = log10(TL_mm),                # log body size
    Body_elong  = TL_mm / BD_mm,               # elongation
    Body_comp   = BD_mm / BW_mm,               # lateral compression
    Pec_AR      = PH_mm^2 / PA_mm2,            # pectoral fin aspect ratio
    Pec_pos     = PFP_mm / BD_mm,              # pectoral fin position
    CP_throttle = (BD_CP_mm)^2 / CPL_mm^2,     # caudal peduncle throttling
    Caud_AR     = CFL_mm^2 / CFA_mm2,          # caudal fin aspect ratio
    # ── Trophic / foraging (4 traits) ─────────────────────────────────
    Mouth_pos   = HMd_mm / BD_mm,              # mouth position (0=ventral,1=dorsal)
    Oral_gape   = pi * (MH_mm * MW_mm) / 4 / BD_mm^2,   # oral gape surface
    Eye_size    = ED_mm^2 / BD_mm^2            # relative eye size
  )

# =============================================================================
# 3. INTRA-SPECIFIC VARIATION OF TRAITS (PRELUDE TO DAY 5)
# =============================================================================
# CV per species per trait — quick look at within-species variability.
trait_cols <- c("Body_size","Body_elong","Body_comp","Pec_AR","Pec_pos",
                "CP_throttle","Caud_AR","Mouth_pos","Oral_gape","Eye_size")

cv_by_sp <- fish_traits_ind |>
  pivot_longer(all_of(trait_cols), names_to = "trait", values_to = "value") |>
  group_by(species, trait) |>
  summarise(mean = mean(value, na.rm = TRUE),
            sd   = sd(value,   na.rm = TRUE),
            # Guard against division by zero / Inf when mean is 0 or near 0
            cv   = ifelse(abs(mean) < 1e-10, NA_real_, sd / abs(mean)),
            .groups = "drop")

p_cv <- ggplot(cv_by_sp, aes(trait, cv, fill = trait)) +
  geom_boxplot(alpha = 0.85, show.legend = FALSE) +
  scale_fill_viridis_d(option = "viridis") +
  labs(title = "Within-species coefficient of variation across the 10 traits",
       x = NULL, y = "CV") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("outputs/day4_within_species_CV.png", p_cv,
       width = 8, height = 4.5, dpi = 300)

# Aggregate to species means (for FD computation downstream)
fish_sp_means_practical <- fish_traits_ind |>
  group_by(species) |>
  summarise(across(all_of(trait_cols), mean, na.rm = TRUE), .groups = "drop")

write_csv(fish_sp_means_practical,
          "outputs/day4_practical_fish_sp_means.csv")

# =============================================================================
# 4. WORK WITH THE FULL 30-SPECIES TRAIT MATRIX
# =============================================================================
# For functional diversity at the community scale we need ALL species in the
# Amazonian dataset. Use the pre-compiled species means (see Day 1/3).

fish_traits_sp <- read_csv("data/fish_traits_sp_means.csv",
                           show_col_types = FALSE) |>
  column_to_rownames("species")

fish_comm <- read_csv("data/fish_communities_amazon.csv",
                      show_col_types = FALSE) |>
  column_to_rownames("site_id")

site_env  <- read_csv("data/fish_site_env.csv", show_col_types = FALSE)
sp_meta   <- read_csv("data/fish_species_metadata.csv",
                      show_col_types = FALSE)

# Align species columns of fish_comm with rows of trait matrix
stopifnot(setequal(colnames(fish_comm), rownames(fish_traits_sp)))
fish_traits_sp <- fish_traits_sp[colnames(fish_comm), , drop = FALSE]

# Trait correlation matrix — useful before FD computation
# corrplot draws on the current graphics device. To save the figure, open
# a png() / pdf() device first then call corrplot.
png("outputs/day4_trait_correlation.png", width = 1500, height = 1500, res = 200)
corrplot::corrplot(
  cor(as.matrix(fish_traits_sp), use = "pairwise.complete.obs"),
  method = "color", type = "upper", addCoef.col = "black",
  number.cex = 0.7, tl.col = "#0A4A5A", tl.srt = 35,
  col = colorRampPalette(c("#2C7BB6","white","#D7191C"))(200))
dev.off()

# Scale traits (different units).
# Guard: drop / warn if any trait has zero variance (would produce NaN columns).
trait_mat_raw <- as.matrix(fish_traits_sp)
trait_sd <- apply(trait_mat_raw, 2, sd, na.rm = TRUE)
if (any(trait_sd == 0 | is.na(trait_sd))) {
  warning("Traits with zero or NA variance dropped before scaling: ",
          paste(names(trait_sd)[trait_sd == 0 | is.na(trait_sd)],
                collapse = ", "))
  trait_mat_raw <- trait_mat_raw[, trait_sd > 0 & !is.na(trait_sd), drop = FALSE]
}
trait_mat_sc <- scale(trait_mat_raw)

# =============================================================================
# 5. COMPUTE FUNCTIONAL DIVERSITY INDICES — dbFD()
# =============================================================================
# Adaptive m: number of PCoA axes can never exceed (n_species - 1).
# Wrapped in tryCatch in case any community has too few species for FRic.
m_use <- min(4, ncol(trait_mat_sc), nrow(trait_mat_sc) - 1)

fd_fish <- tryCatch(
  dbFD(x = trait_mat_sc, a = as.matrix(fish_comm),
       calc.FRic = TRUE, calc.FDiv = TRUE,
       calc.CWM  = TRUE, CWM.type = "dom",
       stand.x   = FALSE, m = m_use, print.pco = FALSE,
       messages  = FALSE),
  error = function(e) {
    message("dbFD() failed with m = ", m_use, " — retrying without FRic.")
    dbFD(x = trait_mat_sc, a = as.matrix(fish_comm),
         calc.FRic = FALSE, calc.FDiv = FALSE,
         calc.CWM  = TRUE,  CWM.type  = "dom",
         stand.x   = FALSE, print.pco = FALSE, messages = FALSE)
  }
)

# Provide NA columns if calc.FRic / calc.FDiv didn't run
if (is.null(fd_fish$FRic)) fd_fish$FRic <- rep(NA_real_, nrow(fish_comm))
if (is.null(fd_fish$FDiv)) fd_fish$FDiv <- rep(NA_real_, nrow(fish_comm))

fd_df <- tibble(site_id = rownames(fish_comm),
                SR    = fd_fish$nbsp,
                FRic  = fd_fish$FRic,
                FEve  = fd_fish$FEve,
                FDiv  = fd_fish$FDiv,
                FDis  = fd_fish$FDis,
                RaoQ  = fd_fish$RaoQ) |>
  left_join(site_env, by = "site_id")

write_csv(fd_df, "outputs/day4_fd_per_site.csv")

# CWM matrix (sites × traits) — left_join may drop duplicate names if a trait
# accidentally shares a name with an env variable; we use suffix() to be safe.
cwm_df <- as.data.frame(fd_fish$CWM) |>
  rownames_to_column("site_id") |>
  left_join(site_env, by = "site_id", suffix = c("_cwm", "_env"))
write_csv(cwm_df, "outputs/day4_cwm_per_site.csv")

# =============================================================================
# 6. FD vs. ENVIRONMENTAL DEGRADATION
# =============================================================================
p_fric <- ggplot(fd_df, aes(degradation, FRic)) +
  geom_smooth(method = "lm", colour = "#0A4A5A", fill = "#A8D5BA") +
  geom_point(aes(colour = river), size = 4) +
  geom_text_repel(aes(label = site_id), size = 3) +
  scale_colour_brewer(palette = "Set2") +
  labs(title = "Functional richness along the Amazonian degradation gradient",
       x = "Degradation index (composite 0–10)",
       y = "FRic")

p_fdis <- ggplot(fd_df, aes(degradation, FDis)) +
  geom_smooth(method = "lm", colour = "#0A4A5A", fill = "#A8D5BA") +
  geom_point(aes(colour = river), size = 4) +
  geom_text_repel(aes(label = site_id), size = 3) +
  scale_colour_brewer(palette = "Set2") +
  labs(title = "Functional dispersion along the Amazonian degradation gradient",
       x = "Degradation index (composite 0–10)", y = "FDis")

p_fd <- p_fric + p_fdis + plot_layout(guides = "collect")
ggsave("outputs/day4_FD_vs_degradation.png", p_fd,
       width = 12, height = 5, dpi = 300)

# =============================================================================
# 7. FUNCTIONAL SPACE — PCoA + species labels coloured by trophic guild
# =============================================================================
sp_dist <- dist(trait_mat_sc)
pcoa    <- ape::pcoa(sp_dist)
# Guard: pcoa() may return < 2 valid axes if data is degenerate
n_axes_avail <- ncol(pcoa$vectors)
n_axes_use   <- min(2, n_axes_avail)
pcoa_df <- as.data.frame(pcoa$vectors[, seq_len(n_axes_use), drop = FALSE]) |>
  rownames_to_column("species")
# Always rename to PC1 (and PC2 if available)
pc_names <- paste0("PC", seq_len(n_axes_use))
colnames(pcoa_df)[2:(1 + n_axes_use)] <- pc_names
if (n_axes_use < 2) pcoa_df$PC2 <- 0   # so plot still works
pcoa_df <- left_join(pcoa_df, sp_meta, by = "species")

p_space <- ggplot(pcoa_df, aes(PC1, PC2, colour = trophic_guild)) +
  geom_hline(yintercept = 0, lty = 2, colour = "grey60") +
  geom_vline(xintercept = 0, lty = 2, colour = "grey60") +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_text_repel(aes(label = species), size = 2.6, max.overlaps = 30) +
  scale_colour_brewer(palette = "Dark2") +
  labs(title = "Amazonian fish functional space (PCoA on 10 morphological traits)",
       subtitle = "Coloured by trophic guild — see Villéger 2017 trait set",
       colour = "Trophic guild")

ggsave("outputs/day4_functional_space.png", p_space,
       width = 10, height = 6.5, dpi = 300)

# =============================================================================
# 8. CWM ALONG ENVIRONMENTAL GRADIENTS
# =============================================================================
cwm_long <- cwm_df |>
  pivot_longer(any_of(colnames(fish_traits_sp)),
               names_to = "trait", values_to = "CWM_value")

p_cwm <- ggplot(cwm_long, aes(degradation, CWM_value)) +
  geom_smooth(method = "lm", colour = "#0A4A5A", fill = "#A8D5BA") +
  geom_point(colour = "#B23A48", alpha = 0.7) +
  facet_wrap(~ trait, scales = "free_y", ncol = 5) +
  labs(title = "CWM of morphological traits along the degradation gradient",
       x = "Degradation index", y = "CWM")
ggsave("outputs/day4_CWM_gradient.png", p_cwm,
       width = 14, height = 6, dpi = 300)

# =============================================================================
# 9. SAVE OBJECTS FOR THE EXERCISES
# =============================================================================
saveRDS(list(fish_traits_ind = fish_traits_ind,
             fish_traits_sp  = fish_traits_sp,
             fish_comm       = fish_comm,
             trait_mat_sc    = trait_mat_sc,
             fd_df           = fd_df,
             cwm_df          = cwm_df,
             pcoa            = pcoa,
             sp_meta         = sp_meta),
        "outputs/day4_objects.rds")

cat("\n=== DAY 4 — outputs written to outputs/ ===\n")
cat("Number of sites:", nrow(fish_comm),
    "  Number of species:", ncol(fish_comm), "\n")
cat("Mean species richness per site:", round(mean(fd_df$SR), 1),
    "  Mean FRic:", round(mean(fd_df$FRic, na.rm=TRUE), 3), "\n")

# =============================================================================
# EXERCISES — see exercises_day4.R
#   Q1 Measurement quality control — TEM, ICC, landmark error
#   Q2 Trait–ecology relationships — Mouth_pos × Body_elong, correlation matrix, PCoA
#   Q3 FD across sites           — degradation gradient, FRic vs FDis sensitivity
#   Q4 From morphology to function — CWM(Oral_gape) vs invertebrates, fish FUSE
# =============================================================================
