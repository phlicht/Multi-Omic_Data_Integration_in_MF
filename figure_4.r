library(pheatmap)
library(ggpubr)
library(ggplotify)
library(dplyr)

# Set Annotation Colors
ann_colors = list(
  Tissue = c(
    "Blood (PBMC)" = "#FFB90F",
    "Nonlesional Skin" = "#0000CD",
    "Patch Skin" = "#8F9C1B",
    "Plaque Skin" = "#2D5A27"
  )
)


###############
## Figure 4A ##
###############

# MetaPhlAn detected 157 viruses
# However, the vast majority is prevalent in only very few samples
# Therefore, we filter the list of viruses to retain the top 20 viruses with the highest variation between samples
top20_var <- apply(metaphlan_virus_inclUnkown, 1, var)
top20_var <- names(sort(top20_var, decreasing = TRUE)[1:20])

superficial_skin_viruses <- superficial_skin_viruses %>%
    dplyr::filter(
      rownames(superficial_skin_viruses) %in% top20_var
    )

# replace 0 with NA
superficial_skin_viruses[superficial_skin_viruses == 0] = NA

# Plot the Heatmap
viruses_heatmap_superficial_skin <- pheatmap(
  superficial_skin_viruses,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = metadata %>% select(Tissue),
  annotation_colors = ann_colors,
  border_color = "NA",
  main = "Superficial Skin"
)


###############
## Figure 4B ##
###############

viruses_heatmap_whole_skin <- pheatmap(
  whole_skin_viruses,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = metadata %>% select(Tissue),
  annotation_colors = ann_colors,
  border_color = "NA",
  main = "Whole Skin"
)

#######################
## Arrange the Plots ##
#######################

figure_4 <- ggarrange(
  ggplotify::as.ggplot(viruses_heatmap_superficial_skin),
  ggplotify::as.ggplot(viruses_heatmap_whole_skin),
  labels = "AUTO"
)


#######################
## Figures 4C and 4D ##
#######################

# Figures 4C and 4D were created using GraphPad Prism.