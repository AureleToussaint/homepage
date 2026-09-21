# ==============================================================
# Master FBE 2026 - Group C - Null models and assembly rules (Oriental + Australian)
# Starter script matching FBE_2026_exercise_C.docx (tested with R 4.3, FD 1.0.12.5)
# Set the working directory to your group folder before running.
# ==============================================================

library(FD); library(vegan); library(ggplot2); library(dplyr)

traits <- read.csv("fish_traits.csv",      stringsAsFactors = FALSE)
comm   <- read.csv("fish_communities.csv", row.names = 1, check.names = FALSE)
env    <- read.csv("env_data.csv",         row.names = 1, stringsAsFactors = FALSE)
rownames(traits) <- traits$species

# Keep only the basins of your two realms and the species present there
my_realms <- c("Oriental", "Australian")
sites     <- rownames(env)[env$Realm %in% my_realms]
comm_sub  <- comm[sites, ]
comm_sub  <- comm_sub[, colSums(comm_sub) > 0]          # drop absent species
tr_sub    <- traits[colnames(comm_sub), 3:9]            # the 7 morphological traits
realm_sp  <- traits[colnames(comm_sub), "Realm"]
cat("Basins:", nrow(comm_sub), " Species:", ncol(comm_sub), "\n")

# ---- Step 1 — Observed functional richness (20 min) ----

set.seed(42)
fd_obs <- dbFD(as.matrix(tr_sub), as.matrix(comm_sub), m = 3,
               calc.FRic = TRUE, calc.CWM = FALSE, print.pco = FALSE)
df <- data.frame(basin = rownames(comm_sub), S = specnumber(comm_sub),
                 FRic = fd_obs$FRic, FDis = fd_obs$FDis,
                 Realm = env[rownames(comm_sub), "Realm"],
                 temp  = env[rownames(comm_sub), "temp_mean"])
print(head(df), digits = 3)

cor.test(df$S, df$FRic, method = "spearman")
cor.test(df$S, df$FDis, method = "spearman")
ggplot(df, aes(S, FRic, colour = Realm)) + geom_point(size = 2.5) +
  geom_smooth(method = "lm", colour = "black", linetype = 2) +
  scale_colour_manual(values = c("#27ae60", "#c0392b")) +
  labs(x = "Species richness S", y = "Observed FRic") + theme_bw()

# ---- Step 2 — Null model and standardised effect sizes (40 min) ----

set.seed(42)
null_arr <- replicate(199, {
  tr_null <- tr_sub[sample(nrow(tr_sub)), , drop = FALSE]
  rownames(tr_null) <- rownames(tr_sub)
  fd_n <- tryCatch(dbFD(as.matrix(tr_null), as.matrix(comm_sub), m = 3,
                        calc.FRic = TRUE, calc.CWM = FALSE,
                        print.pco = FALSE, messages = FALSE),
                   error = function(e) NULL)
  if (is.null(fd_n)) matrix(NA, 2, nrow(comm_sub))
  else rbind(FRic = fd_n$FRic, FDis = fd_n$FDis)
})
null_FRic <- null_arr[1, , ]      # basins x permutations
null_FDis <- null_arr[2, , ]

ses <- function(obs, null_mat)
  (obs - rowMeans(null_mat, na.rm = TRUE)) / apply(null_mat, 1, sd, na.rm = TRUE)
pv  <- function(obs, null_mat) {              # two-sided permutation p-value
  hi <- rowMeans(sweep(null_mat, 1, obs, ">="), na.rm = TRUE)
  lo <- rowMeans(sweep(null_mat, 1, obs, "<="), na.rm = TRUE)
  pmin(1, 2 * pmin(hi, lo)) }
df$SES_FRic <- ses(df$FRic, null_FRic); df$p_FRic <- pv(df$FRic, null_FRic)
df$SES_FDis <- ses(df$FDis, null_FDis); df$p_FDis <- pv(df$FDis, null_FDis)
df$pattern  <- ifelse(df$SES_FRic > 1.96, "overdispersed",
               ifelse(df$SES_FRic < -1.96, "clustered", "random"))
print(df[, c("basin", "Realm", "S", "SES_FRic", "p_FRic", "SES_FDis", "pattern")], digits = 2)

ggplot(df, aes(Realm, SES_FRic, fill = Realm)) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 2, colour = "red") +
  geom_hline(yintercept = 0) +
  geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.08, alpha = 0.5) +
  scale_fill_manual(values = c("#27ae60", "#c0392b")) +
  labs(y = "SES-FRic", title = "Functional assembly patterns") +
  theme_bw() + theme(legend.position = "none")
ggsave("GroupC_SES.png", dpi = 150, width = 6, height = 5)
table(df$Realm, df$pattern)

cor.test(df$temp, df$SES_FRic, method = "spearman")
ggplot(df, aes(temp, SES_FRic, colour = Realm)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 2.5) + geom_smooth(method = "lm", colour = "black") +
  scale_colour_manual(values = c("#27ae60", "#c0392b")) +
  labs(x = "Mean water temperature (°C)", y = "SES-FRic") + theme_bw()
ggsave("GroupC_SES_temp.png", dpi = 150, width = 7, height = 5)

# ---- Step 3 — Robustness: does the conclusion depend on the index? (20 min) ----

cor.test(df$SES_FRic, df$SES_FDis)
df_long <- tidyr::pivot_longer(df, c(SES_FRic, SES_FDis), names_to = "index", values_to = "SES")
ggplot(df_long, aes(Realm, SES, fill = Realm)) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 2, colour = "red") + geom_hline(yintercept = 0) +
  geom_boxplot(alpha = 0.7) + facet_wrap(~ index) +
  scale_fill_manual(values = c("#27ae60", "#c0392b")) +
  theme_bw() + theme(legend.position = "none")
ggsave("GroupC_SES_compare.png", dpi = 150, width = 8, height = 5)
