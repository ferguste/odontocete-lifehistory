### 8. Interpret Results

# Examine the simplex projection of life-history archetypes with family-level overlays
Labels <- c("Bet-hedging", "Reproductive", "Beaked")
simplex <- simplexplot(a3, show_points = FALSE, labels = Labels)

# Plot odontocete families using color overlays
family_colors <- c("Ziphiidae" = "brown", "Physeteridae" = "orange", 
                   "Kogiidae" = "red", "Monodontidae" = "green", 
                   "Delphinidae" = "black", "Phocoenidae" = "blue")

for (fam in names(family_colors)) {
  fam_idx <- which(dat$fam == fam)
  points(simplex$proj_h[fam_idx, ], 
         col = adjustcolor(family_colors[fam], alpha.f = 0.5), 
         pch = 20)
}
legend("topleft", bty = "n", ncol = 2, legend = names(family_colors), 
       fill = family_colors)

# ----- Load and clean phylogenetic tree -----
tree <- read.nexus("C:/Rdata/tree.nex")
tree <- tree[-which(names(tree) == "node.label")]
class(tree) <- "phylo"
plot(tree)

# Confirm tree and dataset match
stopifnot(all(dat$sp %in% tree$tip.label), all(tree$tip.label %in% dat$sp))

# ----- Prepare Data -----
Dataice <- cbind(dat[, c("wtresid", "gestresid", "interresid", "neoresid", 
                         "longresid", "asmresid", "lbathy", "Ice15sMean", "ablat")],
                 lice = log(dat$Ice15sMean + 0.1))
rownames(Dataice) <- dat$sp
summary(Dataice)

# Prepare additional dataset for SEM/correlations
Data2 <- cbind(dat[, c("gestresid", "longresid", "wtresid", "asmresid", 
                       "neoresid", "interresid", "bathymean", "Ice15sMean")],
               latx = abs(dat$latmean) / 90)
rownames(Data2) <- dat$sp

# Rename columns for clarity
Data3 <- Data2 %>% 
  rename(gest = gestresid, long = longresid, wt = wtresid, asm = asmresid,
         neo = neoresid, int = interresid, bathy = bathymean, ice = Ice15sMean, lat = latx)
summary(Data3)

# ----- Correlation Matrix -----
library(GGally)
library(ggplot2)

small <- Dataice %>%
  mutate(latx = abs(ablat) / 90) %>%
  select(wtresid, gestresid, interresid, longresid, neoresid, asmresid, 
         ablat, lbathy, lice)

ggcorr(small, label = TRUE, label_size = 3, color = "grey50")

# ----- Structural Equation Modeling (SEM) -----
library(phylosem)

small <- small %>% 
  rename(gest = gestresid, long = longresid, wt = wtresid, asm = asmresid,
         neo = neoresid, int = interresid, lat = ablat,
         bathy = lbathy, ice = lice)

# Define SEM model
semx <- "
  bathy -> neo, d1
  bathy -> gest, d2
  bathy -> asm, d3
  bathy -> wt, d4
  lat   -> wt,  c1
  lat   -> asm, c2
  ice   -> gest, e1
  ice   -> long, e2
  ice   -> asm,  e3
  gest  -> int,  a1
  gest  -> long, a2
  asm   -> long, a3
  asm   -> gest, a4
  asm   -> int,  a5
  neo   -> int,  a7
"

mx <- phylosem(sem = semx, data = small, tree = tree, data_labels = dat$sp)
summary(mx)
AIC(mx)

# Plot SEM graph
mypp3 <- as_fitted_DAG(mx)
plot(mypp3)

# ----- Generalized Additive Models (GAMs) -----
library(mgcv)

# GAMs vs latitude (ablat)
gam_lat_models <- list(
  wt    = gam(wtresid ~ s(ablat), data = dat, method = "REML"),
  neo   = gam(neoresid ~ s(ablat), data = dat, method = "REML"),
  asm   = gam(asmresid ~ s(ablat), data = dat, method = "REML"),
  inter = gam(interresid ~ s(ablat), data = dat, method = "REML"),
  long  = gam(longresid ~ s(ablat), data = dat, method = "REML")
)

par(mfrow = c(2, 3))
for (name in names(gam_lat_models)) {
  plot(gam_lat_models[[name]], main = paste("Latitude vs", name), shade = TRUE)
}

# GAMs vs ice
gam_ice_models <- list(
  wt    = gam(wtresid ~ s(lice), data = Dataice, method = "REML"),
  neo   = gam(neoresid ~ s(lice), data = Dataice, method = "REML"),
  asm   = gam(asmresid ~ s(lice), data = Dataice, method = "REML"),
  inter = gam(interresid ~ s(lice), data = Dataice, method = "REML"),
  long  = gam(longresid ~ s(lice), data = Dataice, method = "REML")
)

par(mfrow = c(2, 3))
for (name in names(gam_ice_models)) {
  plot(gam_ice_models[[name]], main = paste("Ice vs", name), shade = TRUE)
}

# GAMs vs bathymetry
gam_bathy_models <- list(
  wt    = gam(wtresid ~ s(lbathy), data = dat, method = "REML"),
  neo   = gam(neoresid ~ s(lbathy), data = dat, method = "REML"),
  asm   = gam(asmresid ~ s(lbathy), data = dat, method = "REML"),
  inter = gam(interresid ~ s(lbathy), data = dat, method = "REML"),
  long  = gam(longresid ~ s(lbathy), data = dat, method = "REML")
)

par(mfrow = c(2, 3))
for (name in names(gam_bathy_models)) {
  plot(gam_bathy_models[[name]], main = paste("Bathy vs", name), shade = TRUE)
}

# ----- Phylogenetic Generalized Least Squares (PGLS) -----
library(ape)
library(nlme)
library(phytools)

# Match species between data and tree
dat$species <- as.character(dat$species)
dat <- dat[dat$species %in% tree$tip.label, ]
tree <- drop.tip(tree, setdiff(tree$tip.label, dat$species))
rownames(dat) <- dat$species

# Create phylogenetic correlation structure
corStruct <- corBrownian(phy = tree, form = ~species)

# Create spline basis for latitude
spline_basis <- smoothCon(s(ablat, bs = "tp", k = 10), data = dat)[[1]]
X <- PredictMat(spline_basis, dat)
colnames(X) <- paste0("X", seq_len(ncol(X)))
dat <- cbind(dat, X)

# Run PGLS model using spline basis
gam_latwt_phylo <- gls(wtresid ~ X1 + X2 + X3 + X4 + X5, 
                       data = dat, 
                       correlation = corStruct, 
                       method = "REML")
summary(gam_latwt_phylo)
