# =============================================================================
# DAY 1 — Trait Data: Import, Cleaning and Exploratory Analysis
# Course: Trait-Based Ecology — Theory and Applications in R
# =============================================================================
# OBJECTIVE
#   - Import trait data in long (TRY) and wide format
#   - Detect and handle missing values, outliers and unit issues
#   - Compute species-mean traits
#   - Visualise trait distributions and the Leaf Economics Spectrum
#
# DATA FILES USED
#   data/01b_TRY_format_long.csv  (long format, simulating a TRY export)
#   data/01_species_traits.csv    (wide format, species-mean traits)
#
# PACKAGES REQUIRED
#   tidyverse, ggplot2, GGally, ggrepel, patchwork
# =============================================================================

# ── 0. Setup ──────────────────────────────────────────────────────────────────
library(tidyverse)   # data manipulation + ggplot2
library(GGally)      # pairs plots
library(ggrepel)     # non-overlapping labels
library(patchwork)   # multi-panel figure composition

# Set a consistent ggplot theme for all figures in the course
theme_set(
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "#E1F5EE"),
      strip.text        = element_text(colour = "#085041", face = "bold")
    )
)

# ── 1. Import TRY long-format data ────────────────────────────────────────────
# TRY exports data in long format: one row per observation
# Key columns: AccSpeciesName, TraitID, TraitName, StdValue, UnitName

try_raw <- read_csv("data/01b_TRY_format_long.csv", show_col_types = FALSE)
glimpse(try_raw)

# How many species, traits, and observations?
try_raw |>
  summarise(
    n_obs     = n(),
    n_species = n_distinct(AccSpeciesName),
    n_traits  = n_distinct(TraitID),
    pct_NA    = mean(is.na(StdValue)) * 100
  )

# Which traits are available?
try_raw |>
  count(TraitID, TraitName, UnitName) |>
  print(n = Inf)


# ── 2. Create a concise trait abbreviation ────────────────────────────────────
# TraitName strings from TRY are verbose; we map them to short codes

trait_key <- tribble(
  ~TraitID, ~trait_short, ~unit,
  3117,     "SLA",        "mm2.mg-1",
  47,       "LDMC",       "mg.g-1",
  14,       "LeafN",      "mg.g-1",
  3106,     "Height",     "cm",
  6,        "LSWC",       "g.g-1",
  26,       "SeedMass",   "mg"
)

try_clean <- try_raw |>
  left_join(trait_key, by = "TraitID") |>
  filter(!is.na(StdValue))              # remove missing values (documented)

# How many observations remain?
cat("Observations after removing NA:", nrow(try_clean), "\n")


# ── 3. Outlier detection ───────────────────────────────────────────────────────
# Strategy: flag values beyond mean ± 3 SD *within each trait*
# NOTE: always work on log-scale for right-skewed traits (SLA, SeedMass, Height)

# Visual inspection BEFORE filtering — always look first!
try_clean |>
  ggplot(aes(x = StdValue)) +
  geom_histogram(bins = 40, fill = "#1D9E75", colour = "white", linewidth = 0.3) +
  facet_wrap(~trait_short, scales = "free") +
  labs(title = "Raw trait distributions (before outlier removal)",
       x = "Trait value", y = "Count") +
  theme(axis.text.x = element_text(size = 8))

# Log-transform right-skewed traits for outlier detection
log_traits <- c("SLA", "LDMC", "Height", "SeedMass", "LSWC")

try_filtered <- try_clean |>
  group_by(trait_short) |>
  mutate(
    log_val    = if_else(trait_short %in% log_traits, log(StdValue), StdValue),
    mean_log   = mean(log_val, na.rm = TRUE),
    sd_log     = sd(log_val,   na.rm = TRUE),
    is_outlier = abs(log_val - mean_log) > 3 * sd_log
  ) |>
  filter(!is_outlier) |>
  ungroup()

try_filtered |>
  ggplot(aes(x = log_val)) +
  geom_histogram(bins = 40, fill = "#1D9E75", colour = "white", linewidth = 0.3) +
  facet_wrap(~trait_short, scales = "free") +
  labs(title = "Raw trait distributions (before outlier removal)",
       x = "Trait value", y = "Count") +
  theme(axis.text.x = element_text(size = 8))

# How many observations were removed as outliers?
n_removed <- nrow(try_clean) - nrow(try_filtered)
cat(sprintf("Outliers removed: %d (%.1f%%)\n", n_removed,
            100 * n_removed / nrow(try_clean)))


# ── 4. Compute species-mean traits ────────────────────────────────────────────
# The standard approach: average all valid observations per species × trait
# We also record n_obs (data support) and CV (intraspecific variability signal)

sp_means_long <- try_filtered |>
  group_by(AccSpeciesName, trait_short, unit) |>
  summarise(
    mean_val = mean(StdValue, na.rm = TRUE),
    sd_val   = sd(StdValue,   na.rm = TRUE),
    n_obs    = n(),
    CV       = sd_val / mean_val * 100,   # coefficient of variation (%)
    .groups  = "drop"
  )

# Pivot to wide format: one row per species, one column per trait
sp_traits_wide <- sp_means_long |>
  select(AccSpeciesName, trait_short, mean_val) |>
  pivot_wider(names_from = trait_short, values_from = mean_val)

# How much data coverage do we have?
sp_traits_wide |>
  summarise(across(where(is.numeric), ~mean(!is.na(.)) * 100,
                   .names = "{.col}_pct_complete"))


# ── 5. Load the pre-computed species trait table ───────────────────────────────
# For practicals we also use the clean reference table directly
# (equivalent to the result of the pipeline above)

sp_traits <- read_csv("data/01_species_traits.csv", show_col_types = FALSE)
glimpse(sp_traits)

# Factor levels for PFT (useful for ordered legends later)
sp_traits <- sp_traits |>
  mutate(
    PFT = factor(PFT, levels = c("PFT_stress","PFT_shrub",
                                 "PFT_grass","PFT_geophyte","PFT_ruderal")),
    # Short PFT labels for plotting
    PFT_label = recode(PFT,
      "PFT_stress"   = "Stress-tolerant",
      "PFT_shrub"    = "Dwarf shrub",
      "PFT_grass"    = "Grass",
      "PFT_geophyte" = "Geophyte",
      "PFT_ruderal"  = "Ruderal"
    )
  )


# ── 6. Trait distributions per PFT ────────────────────────────────────────────
# Violin + boxplot combination: best of both worlds
# (violin = density shape, boxplot = median and IQR)

p_sla <- ggplot(sp_traits, aes(x = PFT_label, y = SLA, fill = PFT_label)) +
  geom_violin(alpha = 0.5, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = 21,
               outlier.size = 1.5, colour = "grey30") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = NULL, y = expression(SLA~(mm^2~mg^{-1})),
       title = "A — SLA by plant functional type") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_ldmc <- ggplot(sp_traits, aes(x = PFT_label, y = LDMC, fill = PFT_label)) +
  geom_violin(alpha = 0.5, colour = NA) +
  geom_boxplot(width = 0.18, outlier.shape = 21, colour = "grey30") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = NULL, y = expression(LDMC~(mg~g^{-1})),
       title = "B — LDMC by plant functional type") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_sla + p_ldmc
# INTERPRETATION: Stress-tolerant and shrub species have low SLA and high LDMC
# (conservative strategy); ruderals show the opposite (acquisitive strategy).


# ── 7. The Leaf Economics Spectrum ────────────────────────────────────────────
# Classic LES: SLA positively correlated with Leaf N, negatively with LDMC
# Reference: Wright et al. (2004) Nature

# Panel A: SLA vs. Leaf N
p_les_a <- ggplot(sp_traits, aes(x = log(SLA), y = log(LeafN),
                                 colour = PFT_label, label = species)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey20",
              linewidth = 0.8, linetype = "dashed") +
  ggrepel::geom_text_repel(size = 2.8, max.overlaps = 10,
                           segment.colour = "grey70") +
  scale_colour_brewer(palette = "Set2", name = "PFT") +
  labs(x = "log(SLA)", y = "log(Leaf N)",
       title = "A — LES: SLA vs. Leaf N") +
  theme(legend.position = "bottom")

# Panel B: SLA vs. LDMC (negative arm of LES)
p_les_b <- ggplot(sp_traits, aes(x = log(SLA), y = log(LDMC),
                                 colour = PFT_label)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey20",
              linewidth = 0.8, linetype = "dashed") +
  scale_colour_brewer(palette = "Set2", name = "PFT") +
  labs(x = "log(SLA)", y = "log(LDMC)",
       title = "B — LES: SLA vs. LDMC (negative arm)") +
  theme(legend.position = "bottom")

p_les_a + p_les_b +
  plot_annotation(
    title    = "The Leaf Economics Spectrum — Alpine grassland species",
    subtitle = "30 species across 5 plant functional types",
    caption  = "Each point = one species. Line = OLS regression (log–log scale)."
  )

# Pearson correlation on log-transformed values
with(sp_traits, {
  cat("r(log SLA, log LeafN):", round(cor(log(SLA), log(LeafN)), 3), "\n")
  cat("r(log SLA, log LDMC):", round(cor(log(SLA), log(LDMC)), 3), "\n")
})


# ── 8. Full pairwise trait plot ────────────────────────────────────────────────
# GGally::ggpairs gives a quick multi-trait overview
# Useful to spot unexpected correlations before formal analyses

sp_traits |>
  select(SLA, LDMC, LeafN, Height, LSWC, SeedMass, PFT_label) |>
  mutate(across(where(is.numeric), log)) |>   # log-transform for linearity
  ggpairs(
    aes(colour = PFT_label, alpha = 0.7),
    columns = 1:6,
    upper   = list(continuous = wrap("cor", size = 3)),
    lower   = list(continuous = wrap("points", size = 1.5)),
    diag    = list(continuous = wrap("densityDiag", alpha = 0.4))
  ) +
  labs(title = "Pairwise trait relationships (log-transformed)")


# ── 9. Data quality report ────────────────────────────────────────────────────
# Always document data quality: n observations, CV, data source
# This goes into a supplementary table in a real publication

data_quality <- sp_means_long |>
  group_by(trait_short) |>
  summarise(
    n_species   = n(),
    mean_n_obs  = round(mean(n_obs), 1),
    mean_CV_pct = round(mean(CV, na.rm = TRUE), 1),
    pct_CV_gt30 = round(mean(CV > 30, na.rm = TRUE) * 100, 1)
  )

print(data_quality)
# NOTE: High CV (>30%) for SeedMass is expected — seed mass varies enormously
# among individuals and across years. This is relevant for Day 5 (ITV).


# ── 10. Save clean objects for Day 2 ──────────────────────────────────────────
saveRDS(sp_traits, "outputs/day1_sp_traits_clean.rds")
cat("Saved: outputs/day1_sp_traits_clean.rds\n")

