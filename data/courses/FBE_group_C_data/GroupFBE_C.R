# ══════════════════════════════════════════════════════════
# GROUP C — Null models and assembly rules
# Realms: Oriental + Australian
# Packages: FD, ggplot2, dplyr, vegan
# ══════════════════════════════════════════════════════════
library(FD); library(ggplot2); library(dplyr); library(vegan)

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

my_realms <- c("Oriental", "Australian")
d         <- subset_data(my_realms)
comm_sub  <- d$comm
tr_sub    <- d$traits
cat("Species retained:", ncol(comm_sub), "\nBasins:", nrow(comm_sub), "\n")

# ── Step 1 — Observed FRic ────────────────────────────────
set.seed(42)
fd_obs <- dbFD(as.matrix(tr_sub), as.matrix(comm_sub), m=3,
               calc.FRic=TRUE, calc.CWM=FALSE, print.pco=FALSE)$FRic

df <- data.frame(
  basin=rownames(comm_sub), S=specnumber(comm_sub), FRic=fd_obs,
  Realm=env[rownames(comm_sub),"Realm"],
  temp =env[rownames(comm_sub),"temp_mean"]
)
r_test <- cor.test(df$S, df$FRic, method="spearman")
cat(sprintf("FRic vs S: rho=%.2f p=%.4f\n", r_test$estimate, r_test$p.value))

# ── Step 2 — Null model ───────────────────────────────────
cat("Running null model (199 permutations)...\n")
set.seed(42)
null_FRic <- replicate(199, {
  tr_null <- tr_sub[sample(nrow(tr_sub)), , drop=FALSE]
  rownames(tr_null) <- rownames(tr_sub)
  tryCatch(
    dbFD(as.matrix(tr_null), as.matrix(comm_sub), m=3,
         print.pco=FALSE, calc.CWM=FALSE)$FRic,
    error=function(e) rep(NA_real_, nrow(comm_sub))
  )
}, simplify=TRUE)

null_mean <- apply(null_FRic, 1, mean, na.rm=TRUE)
null_sd   <- apply(null_FRic, 1, sd,   na.rm=TRUE)
SES       <- (fd_obs - null_mean) / null_sd
pval      <- rowMeans(sweep(null_FRic, 1, fd_obs, ">="), na.rm=TRUE)

df$SES     <- round(SES, 2)
df$pval    <- round(pval, 3)
df$pattern <- ifelse(SES > 1.96,"overdispersed",
               ifelse(SES < -1.96,"clustered","random"))

cat("\nResults:\n"); print(df[,c("basin","Realm","S","FRic","SES","pval","pattern")])

ggplot(df, aes(Realm,SES,fill=Realm)) +
  geom_hline(yintercept=c(-1.96,1.96), linetype=2, colour="red", linewidth=0.7) +
  geom_hline(yintercept=0) +
  geom_boxplot(alpha=0.7, outlier.shape=21) +
  geom_jitter(width=0.08, alpha=0.5, size=1.5) +
  scale_fill_manual(values=c("#c0392b","#27ae60")) +
  labs(y="SES-FRic", title="Functional assembly patterns") +
  theme_bw() + theme(legend.position="none")
ggsave("GroupC_SES.png", dpi=150, width=6, height=5)

# ── Step 3 — Summary ─────────────────────────────────────
cat("\nAssembly patterns:\n"); print(table(df$Realm, df$pattern))
r2 <- cor.test(df$temp, df$SES, method="spearman")
cat(sprintf("SES vs temp: rho=%.2f p=%.4f\n", r2$estimate, r2$p.value))
ggplot(df, aes(temp,SES,colour=Realm)) +
  geom_hline(yintercept=0, linetype=2) +
  geom_point(size=2.5) + geom_smooth(method="lm") +
  scale_colour_manual(values=c("#c0392b","#27ae60")) +
  labs(x="Mean temp (°C)", y="SES-FRic") + theme_bw()
ggsave("GroupC_SES_temp.png", dpi=150, width=7, height=5)
