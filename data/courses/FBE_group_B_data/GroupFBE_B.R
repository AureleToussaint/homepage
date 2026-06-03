# ══════════════════════════════════════════════════════════
# GROUP B — FD indices and trait–environment links
# Realms: Palearctic + Ethiopian
# Packages: FD, ggplot2, dplyr, tidyr, vegan
# ══════════════════════════════════════════════════════════
library(FD); library(ggplot2); library(dplyr); library(tidyr); library(vegan)

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

my_realms <- c("Palearctic", "Ethiopian")
d         <- subset_data(my_realms)
comm_sub  <- d$comm
tr_sub    <- d$traits
cat("Species retained:", ncol(comm_sub), "\nBasins:", nrow(comm_sub), "\n")

# ── Step 1 — CWM ─────────────────────────────────────────
cwm       <- functcomp(as.matrix(tr_sub), as.matrix(comm_sub))
cwm$Realm <- env[rownames(cwm), "Realm"]
cwm$temp  <- env[rownames(cwm), "temp_mean"]

ggplot(cwm, aes(temp, BL_mean, colour=Realm)) +
  geom_point(size=2.5, alpha=0.8) +
  geom_smooth(method="lm", colour="black", linetype=2) +
  scale_colour_manual(values=c("#e67e22","#27ae60")) +
  labs(x="Mean temp (°C)", y="CWM body length (cm)",
       title="Bergmann rule? CWM body length vs temperature") + theme_bw()
ggsave("GroupB_CWM.png", dpi=150, width=7, height=5)
m <- lm(BL_mean~temp, data=cwm); print(summary(m))

# ── Step 2 — FD indices ───────────────────────────────────
set.seed(42)
fd <- dbFD(x=as.matrix(tr_sub), a=as.matrix(comm_sub), m=3,
           calc.FRic=TRUE, calc.CWM=FALSE, print.pco=FALSE)

fd_df <- data.frame(
  basin=rownames(comm_sub), S=specnumber(comm_sub),
  FRic=fd$FRic, FEve=fd$FEve, FDiv=fd$FDiv, FDis=fd$FDis, RaoQ=fd$RaoQ,
  Realm=env[rownames(comm_sub),"Realm"],
  temp =env[rownames(comm_sub),"temp_mean"]
)

r_test <- cor.test(fd_df$S, fd_df$FRic, method="spearman")
cat(sprintf("\nFRic vs S: rho=%.2f p=%.4f  (richness confound!)\n",
            r_test$estimate, r_test$p.value))

fd_long <- pivot_longer(fd_df, cols=c(FRic,FEve,FDiv,FDis),
                        names_to="index", values_to="value")
ggplot(fd_long, aes(Realm,value,fill=Realm)) + geom_boxplot(alpha=0.7) +
  facet_wrap(~index, scales="free_y") +
  scale_fill_manual(values=c("#e67e22","#27ae60")) +
  theme_bw() + theme(legend.position="none")
ggsave("GroupB_FDindices.png", dpi=150, width=8, height=6)

# ── Step 3 — FDis vs temperature ─────────────────────────
r2 <- cor.test(fd_df$temp, fd_df$FDis, method="spearman")
cat(sprintf("FDis vs temp: rho=%.2f p=%.4f\n", r2$estimate, r2$p.value))
ggplot(fd_df, aes(temp,FDis,colour=Realm)) +
  geom_point(size=2.5) + geom_smooth(method="lm") +
  scale_colour_manual(values=c("#e67e22","#27ae60")) +
  labs(x="Mean temp (°C)", y="FDis") + theme_bw()
ggsave("GroupB_FDis_temp.png", dpi=150, width=7, height=5)
