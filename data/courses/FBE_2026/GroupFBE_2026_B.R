# ==============================================================
# Master FBE 2026 - Group B - FD indices and trait-environment links (Palearctic + Ethiopian)
# Starter script matching FBE_2026_exercise_B.docx (tested with R 4.3, FD 1.0.12.5)
# Set the working directory to your group folder before running.
# ==============================================================

library(FD); library(vegan); library(ggplot2); library(dplyr); library(tidyr)

traits <- read.csv("fish_traits.csv",      stringsAsFactors = FALSE)
comm   <- read.csv("fish_communities.csv", row.names = 1, check.names = FALSE)
env    <- read.csv("env_data.csv",         row.names = 1, stringsAsFactors = FALSE)
rownames(traits) <- traits$species

# Keep only the basins of your two realms and the species present there
my_realms <- c("Palearctic", "Ethiopian")
sites     <- rownames(env)[env$Realm %in% my_realms]
comm_sub  <- comm[sites, ]
comm_sub  <- comm_sub[, colSums(comm_sub) > 0]          # drop absent species
tr_sub    <- traits[colnames(comm_sub), 3:9]            # the 7 morphological traits
realm_sp  <- traits[colnames(comm_sub), "Realm"]
cat("Basins:", nrow(comm_sub), " Species:", ncol(comm_sub), "\n")

# ---- Step 1 — Community-weighted means (20 min) ----

cwm       <- functcomp(as.matrix(tr_sub), as.matrix(comm_sub))
cwm$Realm <- env[rownames(cwm), "Realm"]
cwm$temp  <- env[rownames(cwm), "temp_mean"]
print(head(cwm), digits = 3)
cwm |> group_by(Realm) |> summarise(across(BL_mean:CaudalPed_ratio, mean))

m1 <- lm(BL_mean ~ temp, data = cwm); summary(m1)
ggplot(cwm, aes(temp, BL_mean, colour = Realm)) +
  geom_point(size = 2.5) + geom_smooth(method = "lm", colour = "black", linetype = 2) +
  scale_colour_manual(values = c("#e67e22", "#27ae60")) +
  labs(x = "Mean water temperature (°C)", y = "CWM body length (cm)") + theme_bw()
ggsave("GroupB_CWM.png", dpi = 150, width = 7, height = 5)
# Same test for a trophic trait
summary(lm(MouthPos ~ temp, data = cwm))

# ---- Step 2 — Functional diversity indices (40 min) ----

set.seed(42)
fd <- dbFD(x = as.matrix(tr_sub), a = as.matrix(comm_sub), m = 3,
           calc.FRic = TRUE, calc.CWM = FALSE, print.pco = FALSE)
fd_df <- data.frame(basin = rownames(comm_sub), S = specnumber(comm_sub),
                    FRic = fd$FRic, FEve = fd$FEve, FDiv = fd$FDiv,
                    FDis = fd$FDis, RaoQ = fd$RaoQ,
                    Realm = env[rownames(comm_sub), "Realm"],
                    temp  = env[rownames(comm_sub), "temp_mean"])
print(head(fd_df), digits = 3)

idx <- c("FRic", "FEve", "FDiv", "FDis", "RaoQ")
sapply(idx, function(i) {
  ct <- cor.test(fd_df$S, fd_df[[i]], method = "spearman")
  c(rho = round(ct$estimate, 2), p = round(ct$p.value, 4)) })
ggplot(fd_df, aes(S, FRic, colour = Realm)) + geom_point(size = 2.5) +
  geom_smooth(method = "lm", colour = "black", linetype = 2) +
  scale_colour_manual(values = c("#e67e22", "#27ae60")) +
  labs(x = "Species richness S", y = "FRic") + theme_bw()

fd_long <- pivot_longer(fd_df, cols = all_of(idx), names_to = "index", values_to = "value")
ggplot(fd_long, aes(Realm, value, fill = Realm)) +
  geom_boxplot(alpha = 0.7) + facet_wrap(~ index, scales = "free_y") +
  scale_fill_manual(values = c("#e67e22", "#27ae60")) +
  theme_bw() + theme(legend.position = "none")
ggsave("GroupB_FDindices.png", dpi = 150, width = 8, height = 6)
sapply(idx, function(i) round(wilcox.test(fd_df[[i]] ~ fd_df$Realm)$p.value, 4))

# ---- Step 3 — Trait–environment relationships (20 min) ----

fd_df <- cbind(fd_df, env[rownames(comm_sub), c("discharge", "altitude", "pH")])
cor.test(fd_df$temp, fd_df$FDis, method = "spearman")
m2 <- lm(FDis ~ temp + log(discharge) + altitude, data = fd_df); summary(m2)
ggplot(fd_df, aes(temp, FDis, colour = Realm)) + geom_point(size = 2.5) +
  geom_smooth(method = "lm", colour = "black") +
  scale_colour_manual(values = c("#e67e22", "#27ae60")) +
  labs(x = "Mean water temperature (°C)", y = "FDis") + theme_bw()
ggsave("GroupB_FDis_temp.png", dpi = 150, width = 7, height = 5)
