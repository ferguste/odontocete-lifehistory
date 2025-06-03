## ---------------------------------------------
## SECTION 1: Load Required Libraries
## ---------------------------------------------
library(dplyr)
library(ggplot2)
library(ggdendro)
library(mvMORPH)
library(pglm)
library(MCMCglmm)
library(phytools)
library(ape)
library(nlme)
library(archetypes)

## ---------------------------------------------
## SECTION 2: Prepare Data for Clustering
## ---------------------------------------------
# Select residuals for clustering
clus_W <- dat %>%
  select(sp, wtresid, neoresid, longresid, asmresid, interresid, gestresid)

# Standardize the residuals (robust to outliers)
wcars.use <- clus_W[, -1]
wmedians <- apply(wcars.use, 2, median)
wmads <- apply(wcars.use, 2, mad)
wcars.use <- scale(wcars.use, center = wmedians, scale = wmads)
rownames(wcars.use) <- clus_W$sp

# Calculate distance and perform hierarchical clustering
wcars.dist <- dist(wcars.use)
cars.hclustW <- hclust(wcars.dist)

## ---------------------------------------------
## SECTION 3: Plot Dendrogram for 2 Clusters
## ---------------------------------------------
whale_clust <- cutree(cars.hclustW, k = 2)
dendrW <- dendro_data(cars.hclustW, type = "rectangle")
wclust.df <- data.frame(label = names(whale_clust), cluster = factor(whale_clust))
dendrW[["labels"]] <- merge(dendrW[["labels"]], wclust.df, by = "label")

pal <- c('#C70039', '#04294e')

ggplot() + 
  geom_segment(data = segment(dendrW), aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_text(data = label(dendrW), aes(x, y, label = label, hjust = 0, color = cluster), 
            size = 4, family = 'Roboto Medium', fontface = 'italic') +
  scale_color_manual(values = pal) +
  coord_flip() + 
  scale_y_reverse(expand = c(0.5, 0), breaks = rev(seq(0, 12, by = 2))) + 
  labs(y = 'height') +
  theme_minimal(base_family = "Roboto Medium") +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = 'none'
  )

## ---------------------------------------------
## SECTION 4: Plot Dendrogram for 4 Clusters
## ---------------------------------------------
whale_clust <- cutree(cars.hclustW, k = 4)
dendrW <- dendro_data(cars.hclustW, type = "rectangle")
wclust.df <- data.frame(label = names(whale_clust), cluster = factor(whale_clust))
dendrW[["labels"]] <- merge(dendrW[["labels"]], wclust.df, by = "label")

pal <- c('#C70039', '#04294e', "#ffa500", "#4caf50")

ggplot() + 
  geom_segment(data = segment(dendrW), aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_text(data = label(dendrW), aes(x, y, label = label, hjust = 0, color = cluster), 
            size = 4, family = 'Roboto Medium', fontface = 'italic') +
  scale_color_manual(values = pal) +
  coord_flip() + 
  scale_y_reverse(expand = c(0.5, 0), breaks = rev(seq(0, 12, by = 2))) + 
  labs(y = 'height') +
  theme_minimal(base_family = "Roboto Medium") +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = 'none'
  )

## ---------------------------------------------
## SECTION 5: Phylogenetic ANOVA for Env Variables
## ---------------------------------------------
tree <- drop.tip(tree, setdiff(tree$tip.label, dax$sp))
rownames(dax) <- dax$sp

# Example Phylogenetic ANOVAs
phylANOVA(tree, x = dax$Archetype, y = dax$ablat, nsim = 1000)
phylANOVA(tree, x = dax$Archetype, y = dax$Ice15sMean, nsim = 1000)
phylANOVA(tree, x = dax$Archetype, y = dax$lbathy, nsim = 1000)

## ---------------------------------------------
## SECTION 6: PGLS Models for Environmental Traits
## ---------------------------------------------
corPag <- corPagel(1, phy = tree, fixed = FALSE)

pgls_ablat <- gls(ablat ~ Archetype, data = dax, correlation = corPag, method = "REML")
summary(pgls_ablat)

pgls_ice <- gls(Ice15sMean ~ Archetype, data = dax, correlation = corPag, method = "REML")
summary(pgls_ice)

pgls_bathy <- gls(lbathy ~ Archetype, data = dax, correlation = corPag, method = "REML")
summary(pgls_bathy)

## ---------------------------------------------
## SECTION 7: Archetypes Analysis
## ---------------------------------------------
data <- dat %>%
  select(wtresid, gestresid, longresid, neoresid, asmresid, interresid) %>%
  na.omit()

# Fit archetype models from 1 to 10 archetypes
set.seed(1986)
archetype_model <- stepArchetypes(data, k = 1:10, verbose = TRUE)

# Scree plot to determine optimal number of archetypes
screeplot(archetype_model)

# Choose the 3-archetype model
a3 <- bestModel(archetype_model[[3]])

# Visualization and summaries
summary(a3)
barplot(a3, data)
panorama(a3, data)
pcplot(a3, data)
simplexplot(a3, show_direction = TRUE)

# Assign individuals to archetypes
alphas <- a3$alphas
assigned_archetypes <- apply(alphas, 1, which.max)
table(assigned_archetypes)

# Combine with original dataset
original_data_with_archetypes <- cbind(dat, Archetype = assigned_archetypes)
write.csv(original_data_with_archetypes, "data_with_archetypes.csv", row.names = FALSE)

# Inspect archetype parameters
coef(a3)
kappa(a3, data)
nparameters(a3, data)
parameters(a3, data)
residuals(a3, data)
