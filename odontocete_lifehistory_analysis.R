## odontocete_lifehistory_analysis.R
## Author: Dr. Steven Ferguson, Arctic Region, DFO
## Project Start: 1 March 2023 – Ongoing
## Purpose: Analyze odontocete life history traits using phylogenetic methods

# 1) Load required packages --------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(GGally)
library(ggdendro)
library(dendextend)
library(extrafont)
library(ape)
library(geiger)
library(phytools)
library(nlme)

# 2) Load data ---------------------------------------------
data_path <- "data/"  # relative path for reproducibility

dat <- read.csv(paste0(data_path, "odontocete.csv"))
dax <- read.csv(paste0(data_path, "data_with_archetypes.csv"))
tree <- read.nexus(paste0(data_path, "tree.nex"))

summary(dax)
rownames(dax) <- dax$sp  # species column

# 3) Color dendrogram by Archetype -------------------------
# Convert tree to dendrogram
dend <- as.dendrogram(tree)
label_order <- labels(dend)
dax$Archetype <- factor(dax$Archetype, levels = c(1, 2, 3), 
                        labels = c("Reproducers", "Bet-hedgers", "Beaked whales"))

# Assign colors
group_colors <- c("red", "blue", "forestgreen")
names(group_colors) <- levels(dax$Archetype)
label_archetype <- dax[label_order, "Archetype"]
label_colors <- group_colors[as.character(label_archetype)]

# Color labels
dend_colored <- dend %>%
  set("labels_col", label_colors) %>%
  set("labels_cex", 0.8)

# Convert to ggplot-compatible format
ddata <- dendro_data(dend_colored, type = "rectangle")

# Create label-color lookup
label_df <- data.frame(
  label = ddata$labels$label,
  Archetype = dax[ddata$labels$label, "Archetype"]
)
label_df$color <- group_colors[as.character(label_df$Archetype)]
ddata$labels <- merge(ddata$labels, label_df, by = "label", sort = FALSE)

# Plot dendrogram
dendro_plot <- ggplot() +
  geom_segment(data = ddata$segments, 
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = ddata$labels, 
            aes(x = x, y = y - 0.02, label = label, color = Archetype), 
            hjust = 1, angle = 0, size = 3) +
  scale_color_manual(values = group_colors) +
  labs(y = "Time", x = NULL, title = "Dendrogram with Species Colored by Archetype") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 10),
        legend.position = "right") +
  coord_flip()
