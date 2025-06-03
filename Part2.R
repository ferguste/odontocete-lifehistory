## --------------------------------------------
## Ancestral Trait Reconstructions
## --------------------------------------------
library(phytools)
library(ggplot2)

# Load your tree and trait data
# tree <- read.tree("path/to/tree_file")
# dat <- read.csv("path/to/data.csv")

## Ensure species names match between tree and data
rownames(dat) <- dat$sp
name.check(tree, dat)

## Helper function to reconstruct and plot
plot_ancestral_trait <- function(trait_col, title_text, fig_path) {
  mass <- name_vec(trait_col)
  obj <- contMap(tree, mass, plot = FALSE)
  obj <- setMap(obj, invert = TRUE)
  obj <- setMap(obj, colors = c("yellow", "darkcyan", "purple"))
  
  pdf(fig_path, width = 11, height = 8.5)
  plot.contMap(
    obj, type = "phylogram",
    legend = 0.6 * max(nodeHeights(tree)),
    fsize = c(0.5, 0.7),
    outline = FALSE, lwd = 3, mar = c(2, 2, 2, 2)
  )
  title(title_text, line = 1, cex.main = 1.5)
  add.scale.bar(lwd = 2, cex = 0.7, font = 2, col = "red")
  dev.off()
}

## Fig 1–4: Individual reconstructions
plot_ancestral_trait("ablat", "Fig 1: Latitude", "../figures/fig1_latitude.pdf")
plot_ancestral_trait("lbathy", "Fig 2: Bathymetry", "../figures/fig2_bathymetry.pdf")
plot_ancestral_trait("Ice15sMean", "Fig 3: Sea Ice", "../figures/fig3_ice.pdf")
plot_ancestral_trait("wtresid", "Fig 4: Shape Residual", "../figures/fig4_shape.pdf")

## --------------------------------------------
## Composite Figure: Ancestral Reconstructions
## --------------------------------------------

png("../figures/Figure3_AncestralReconstruction.png", width = 1800, height = 2400, res = 300)
par(mfrow = c(3, 1), mar = c(2, 2, 4, 2))

for (trait_info in list(
  list(col = "ablat", title = "A. Latitude"),
  list(col = "lbathy", title = "B. Bathymetry"),
  list(col = "Ice15sMean", title = "C. Sea Ice")
)) {
  mass <- name_vec(trait_info$col)
  obj <- contMap(tree, mass, plot = FALSE)
  obj <- setMap(obj, invert = TRUE)
  obj <- setMap(obj, colors = c("yellow", "darkcyan", "purple"))
  
  plot.contMap(
    obj, type = "phylogram",
    legend = 0.6 * max(nodeHeights(tree)),
    fsize = c(0.5, 0.7),
    outline = FALSE, lwd = 3, mar = c(2, 2, 2, 2)
  )
  title(trait_info$title, line = 1, cex.main = 1.5)
}
dev.off()

## --------------------------------------------
## Fig 5: Gestation vs. Neonate Size
## --------------------------------------------

p <- ggplot(dat, aes(x = loggest, y = logneo)) +
  geom_point(aes(color = fam), size = 4) +
  stat_ellipse(aes(color = fam), size = 1) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#000000") +
  geom_text(aes(label = sort), size = 3, hjust = -0.5) +
  labs(
    x = expression(log[10]~Gestation),
    y = expression(log[10]~Neonate),
    color = ""
  ) +
  theme_classic(base_size = 16) +
  theme(text = element_text(family = "Roboto Medium"))

ggsave("../figures/fig5_gestation_vs_neonate.pdf", plot = p, height = 8.5, width = 11, device = cairo_pdf)
