###############
## Figure 3A ##
###############
library(MOFA2)
library(ggplot2)

var_explained_per_factor <- plot_variance_explained(
    MOFA2_model.trained,
    factors = 1:5,
    max_r2 = 30
) +
    theme(
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = .9)
    )

###############
## Figure 3B ##
###############
library(MOFA2)
library(pheatmap)

# Italicize the names of the bugs
BugNames_italicized <- MOFA2_model.trained@data$Microbiome$group1 %>% rownames
BugNames_italicized <- lapply(
  BugNames_italicized,
  function(x) bquote(italic(.(x)))
)

# Plot the heatmap
microbiome_weigths_heatmap <- plot_weights_heatmap(
    MOFA2_model.trained,
    view = "Microbiome",
    factors = 1:5,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    labels_col = as.expression(BugNames_italicized)
)

####################
## Figure 3C - 3E ##
####################

# Figure 3C
factor_1_gsea <- plot_enrichment(
    MOFA2_model.trained,
    factor = 1,
    max.pathways = 11
) +
    theme_bw() +
    labs(title = "Factor 1") +
    theme(
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.text = element_text(size = 13),
        title = element_text(size = 13),
        text = element_text(color = "black")
    )

# Figure 3D
factor_4_gsea <- plot_enrichment(
    MOFA2_model.trained,
    factor = 4,
    max.pathways = 9
) +
    theme_bw() +
    labs(title = "Factor 4") +
    theme(
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.text = element_text(size = 13),
        title = element_text(size = 13),
        text = element_text(color = "black")
    )

# Figure 3E
factor_5_gsea <- plot_enrichment(
    MOFA2_model.trained,
    factor = 5,
    max.pathways = 8
) +
    theme_bw() +
    labs(title = "Factor 5") +
    theme(
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.text = element_text(size = 13),
        title = element_text(size = 13),
        text = element_text(color = "black")
    )


#######################
## Arrange the plots ##
#######################
library(ggplot2)
library(ggpubr)
library(ggplotify)

# Arrange Plots 3A and 3B
## Plots 3A and 3B were arranged manually using Microsoft PowerPoint

# Arrange Plots 3C, 3D, and 3E
plot3C_3D <- ggarrange(
    factor_1_gsea,
    factor_4_gsea,
    labels = c("C", "D")
)
plot_3E <- ggarrange(
    factor_5_gsea,
    labels = "E"
)

plot3C_3D_3E <- ggarrange(
    plot3C_3D,
    plot_3E,
    ncol = 1,
    nrow = 2
)

## Plots 3A, 3B, and plot3C_3D_3E were subsequently arranged to figure 3 in Microsoft PowerPoint