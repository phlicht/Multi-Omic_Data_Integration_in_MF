library(MOFA2)
library(ggplot2)
library(dplyr)

# list for saving the output
suppl_material_5_list <- list()

########################
## GSEA Visualization ##
########################

# Loop over the factors
for (n in 1:5) {
    factor <- paste0("Factor ", n)
    print(paste0("Generating GSEA Output for ", factor))

    suppl_material_5_list[["GSEA"]][[factor]] <- plot_enrichment(
        MOFA2_model.trained,
        factor = n,
        max.pathways = 12
    ) +
        theme_bw() +
        labs(title = factor) +
        theme(
            axis.title.y = element_blank(),
            axis.title.x = element_text(size = 13),
            axis.text = element_text(size = 13),
            title = element_text(size = 13),
            text = element_text(color = "black")
        )
}

###################
## GSEA Detailed ##
###################

# Set ggrepel to Inf to enable naming of all selected genes
options(ggrepel.max.overlaps = Inf)

# Loop over the factors
for (n in 1:5) {
    factor <- paste0("Factor ", n)
    print(paste0("Generating GSEA Output for ", factor))

    title_text <- paste0("Top Genes of Enriched Pathways in Factor ", n)

    suppl_material_5_list[["GSEA_DETAILED"]][[factor]] <- plot_enrichment_detailed(
        MOFA2_model.trained,
        factor = n,
        max.genes = 20,
        max.pathways = 8
    ) +
        labs(title = title_text) +
        theme(
            axis.text = element_text(size = 11, color = "black"),
            title = element_text(size = 12),
            text = element_text(color = "black")
        )
}