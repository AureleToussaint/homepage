# ==============================================================
# Master FBE 2026 - Group A - Building and reading a functional space (Neotropical + Nearctic)
# Starter script matching FBE_2026_exercise_A.docx (tested with R 4.3, FD 1.0.12.5)
# Set the working directory to your group folder before running.
# ==============================================================


library(FD); library(vegan); library(ggplot2); library(dplyr); library(reshape2)

traits <- read.csv("fish_traits.csv", stringsAsFactors = FALSE)
rownames(traits) <- traits$species

# Species pool = the species native to your two realms (column Realm of the trait table)
my_realms <- c("Neotropical", "Nearctic")
tr_sub    <- traits[traits$Realm %in% my_realms, 3:9]      # the 7 morphological traits
realm_sp  <- traits$Realm[traits$Realm %in% my_realms]
table(realm_sp)

# ---- Step 1 — Explore the traits (20 min) ----

round(colMeans(!is.na(tr_sub)) * 100, 1)           # % non-NA per trait
ok_sp    <- rowMeans(!is.na(tr_sub)) >= 0.70
tr_clean <- tr_sub[ok_sp, ]
realm_vec <- realm_sp[ok_sp]
summary(tr_clean)

cor_mat  <- cor(tr_clean, use = "pairwise.complete.obs")
cor_melt <- melt(round(cor_mat, 2))
ggplot(cor_melt, aes(Var1, Var2, fill = value, label = value)) +
  geom_tile(colour = "white") + geom_text(size = 2.8) +
  scale_fill_gradient2(low = "#c0392b", mid = "white", high = "#8e44ad", midpoint = 0) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
which(abs(cor_mat) > 0.7 & upper.tri(cor_mat), arr.ind = TRUE)

# ---- Step 2 — Build and read the functional space (40 min) ----

pca     <- prcomp(tr_clean, scale. = TRUE, center = TRUE)
var_pct <- round(100 * summary(pca)$importance[2, ], 1)
cumsum(var_pct)                                       # how many axes for > 70 % ?
round(pca$rotation[, 1:3], 2)                         # loadings
sort(abs(pca$rotation[, 1]), decreasing = TRUE)       # traits driving PC1

scores       <- as.data.frame(pca$x[, 1:2])
scores$Realm <- realm_vec
scores$species <- rownames(tr_clean)
ggplot(scores, aes(PC1, PC2, colour = Realm, fill = Realm)) +
  geom_point(alpha = 0.7, size = 2.5) +
  stat_ellipse(geom = "polygon", alpha = 0.08, linewidth = 0.7) +
  scale_colour_manual(values = c("#2980b9", "#8e44ad")) +
  scale_fill_manual(values   = c("#2980b9", "#8e44ad")) +
  labs(x = paste0("PC1 (", var_pct[1], "%)"), y = paste0("PC2 (", var_pct[2], "%)"),
       title = "Functional space - Neotropical vs Nearctic") + theme_bw()
ggsave("GroupA_funspace.png", dpi = 150, width = 7, height = 5)

realm_mat <- t(sapply(my_realms, function(r) as.numeric(realm_vec == r)))
colnames(realm_mat) <- rownames(tr_clean)
fd_realm <- dbFD(tr_clean, realm_mat, m = 2, stand.FRic = TRUE,
                 calc.FRic = TRUE, calc.CWM = FALSE, print.pco = FALSE)
round(fd_realm$FRic, 3)                                # FRic per realm (0-1)

# Overlap: share of species of one realm lying inside the other realm's convex hull
library(geometry)
in_hull <- function(pts, hull_pts) {
  h <- convhulln(hull_pts); inhulln(h, as.matrix(pts)) }
neo <- as.matrix(scores[scores$Realm == "Neotropical", 1:2])
nea <- as.matrix(scores[scores$Realm == "Nearctic",    1:2])
ov_nea <- 100 * mean(in_hull(nea, neo))    # % Nearctic species inside Neotropical hull
ov_neo <- 100 * mean(in_hull(neo, nea))    # % Neotropical species inside Nearctic hull
round(c(Nearctic_in_Neotropical = ov_nea, Neotropical_in_Nearctic = ov_neo), 1)

# ---- Step 3 — Functional uniqueness (20 min) ----

pc12  <- as.matrix(pca$x[, 1:2])
d_mat <- as.matrix(dist(pc12)); diag(d_mat) <- NA
scores$uniq <- apply(d_mat, 1, min, na.rm = TRUE)
head(scores[order(-scores$uniq), c("species", "Realm", "uniq")], 5)   # 5 most unique

wilcox.test(uniq ~ Realm, data = scores)
ggplot(scores, aes(Realm, uniq, fill = Realm)) +
  geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = c("#2980b9", "#8e44ad")) +
  labs(y = "Functional uniqueness (nearest-neighbour distance)") +
  theme_bw() + theme(legend.position = "none")
ggsave("GroupA_uniqueness.png", dpi = 150, width = 6, height = 5)
