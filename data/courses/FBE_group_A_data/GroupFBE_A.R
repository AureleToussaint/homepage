# ══════════════════════════════════════════════════════════
# GROUP A — Functional space construction + visualisation
# Realms: Neotropical + Nearctic
# Packages: FD, ggplot2, dplyr, reshape2, vegan
# ══════════════════════════════════════════════════════════
library(FD); library(ggplot2); library(dplyr); library(reshape2); library(vegan)

# ── Load data ─────────────────────────────────────────────
traits <- read.csv("fish_traits.csv",       stringsAsFactors = FALSE)
comm   <- read.csv("fish_communities.csv",  row.names = 1, check.names = FALSE)
env    <- read.csv("env_data.csv",          row.names = 1, stringsAsFactors = FALSE)
rownames(traits) <- traits$species

# Helper: subset community + traits, drop species absent from all selected sites
subset_data <- function(my_realms) {
  sites     <- rownames(env[env$Realm %in% my_realms, ])
  sp_all    <- intersect(colnames(comm[sites, ]), traits$species)
  comm_sub  <- comm[sites, sp_all, drop = FALSE]
  comm_sub  <- comm_sub[, colSums(comm_sub) > 0, drop = FALSE]   # remove absent species
  sp_ok     <- colnames(comm_sub)
  tr_sub    <- traits[traits$species %in% sp_ok, 3:9, drop = FALSE]
  list(comm = comm_sub, traits = tr_sub, sites = sites,
       realm_sp = traits$Realm[traits$species %in% sp_ok])
}

my_realms <- c("Neotropical", "Nearctic")
d         <- subset_data(my_realms)
comm_sub  <- d$comm
tr_sub    <- d$traits
realm_sp  <- d$realm_sp
cat("Species retained:", ncol(comm_sub), "\nBasins:", nrow(comm_sub), "\n")

# ── Step 1 — Trait completeness ───────────────────────────
cat("\nCompleteness per trait (%):\n")
print(round(colMeans(!is.na(tr_sub)) * 100, 1))

# Drop species with too many NAs (>30% missing)
ok_sp    <- rowMeans(!is.na(tr_sub)) >= 0.70
tr_clean <- tr_sub[ok_sp, ]
realm_vec <- realm_sp[ok_sp]

# Correlation heatmap
cor_mat  <- cor(tr_clean, use = "pairwise.complete.obs")
cor_melt <- reshape2::melt(round(cor_mat, 2))
ggplot(cor_melt, aes(Var1, Var2, fill = value, label = value)) +
  geom_tile(colour = "white") + geom_text(size = 2.8) +
  scale_fill_gradient2(low="#c0392b",mid="white",high="#8e44ad",midpoint=0) +
  theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1))

# ── Step 2 — PCA + functional space ─────────────────────
pca     <- prcomp(tr_clean, scale.=TRUE, center=TRUE)
var_pct <- round(100 * summary(pca)$importance[2, 1:4], 1)
cat("\nVariance explained: PC1=",var_pct[1],"% PC2=",var_pct[2],"%\n")
cat("PC1 loadings:\n"); print(sort(abs(pca$rotation[,1]),decreasing=TRUE))

scores       <- as.data.frame(pca$x[,1:2])
scores$Realm <- realm_vec

ggplot(scores[scores$Realm %in% my_realms, ],
       aes(PC1, PC2, colour=Realm, fill=Realm)) +
  geom_point(alpha=0.6, size=2.5) +
  stat_ellipse(geom="polygon", alpha=0.08, linewidth=0.7) +
  scale_colour_manual(values=c("#8e44ad","#2980b9")) +
  scale_fill_manual(  values=c("#8e44ad","#2980b9")) +
  labs(title="Functional space — Neotropical vs Nearctic",
       x=paste0("PC1 (",var_pct[1],"%)"), y=paste0("PC2 (",var_pct[2],"%)")) +
  theme_bw()
ggsave("GroupA_funspace.png", dpi=150, width=7, height=5)

# Overlap (% points inside other realm bounding box)
neo_pc <- as.matrix(scores[scores$Realm=="Neotropical",1:2])
nea_pc <- as.matrix(scores[scores$Realm=="Nearctic",   1:2])
bb_overlap <- function(a, b)
  mean(a[,1]>=min(b[,1]) & a[,1]<=max(b[,1]) & a[,2]>=min(b[,2]) & a[,2]<=max(b[,2]))
cat(sprintf("\nOverlap Neotropical-in-Nearctic: %.1f%%\n", 100*bb_overlap(neo_pc,nea_pc)))
cat(sprintf("Overlap Nearctic-in-Neotropical: %.1f%%\n",   100*bb_overlap(nea_pc,neo_pc)))

# ── Step 3 — Functional uniqueness ───────────────────────
pc12  <- as.matrix(pca$x[,1:2])
d_mat <- as.matrix(dist(pc12)); diag(d_mat) <- NA
uniq  <- apply(d_mat, 1, min, na.rm=TRUE)

uniq_df <- data.frame(species=rownames(tr_clean), Realm=realm_vec, uniq=uniq)
cat("\n5 most unique species:\n")
print(head(uniq_df[order(-uniq_df$uniq), c("species","Realm","uniq")], 5))

sub_uniq <- uniq_df[uniq_df$Realm %in% my_realms, ]
wt <- wilcox.test(uniq~Realm, data=sub_uniq)
cat(sprintf("Wilcoxon p=%.4f\n", wt$p.value))
ggplot(sub_uniq, aes(Realm,uniq,fill=Realm)) +
  geom_boxplot(alpha=0.7) + geom_jitter(width=0.1,alpha=0.4,size=1.5) +
  scale_fill_manual(values=c("#8e44ad","#2980b9")) +
  labs(y="Functional uniqueness") + theme_bw() + theme(legend.position="none")
ggsave("GroupA_uniqueness.png", dpi=150, width=6, height=5)
