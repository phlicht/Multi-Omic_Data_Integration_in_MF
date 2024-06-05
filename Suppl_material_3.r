library(DESeq2)
library(pathfindR)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

###################
## Preprocessing ##
###################

# Load Data into DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = metadata,
    design = ~Patient+Stage
)

# Run DESeq2
dds <- DESeq2(dds)

# Extract DEG, shrinked by apeglm
res_Patch <- lfcShrink(dds, coef = "stage_Patch_vs_nonlesional", type = "ageglm")
res_Plaque <- lfcShrink(dds, coef = "stage_Plaque_vs_nonlesional", type = "ageglm")

## Gather informations for pathfindR
DESeq2_DEG_Patch <- res_Patch %>%
    tibble::rownames_to_column("Gene_symbol") %>%
    rename(
        logFC = lfcSE,
        p_val = pvalue
    ) %>%
    select(Gene_symbol, logFC, p_val)

DESeq2_DEG_Plaque <- res_Plaque %>%
    tibble::rownames_to_column("Gene_symbol") %>%
    rename(
        logFC = lfcSE,
        p_val = pvalue
    ) %>%
    select(Gene_symbol, logFC, p_val)


# Transform and extract DESeq2 normalized count ata using variance stabilizing transformation
vsd <- vst(dds)
vsd <- assay(vsd)

# Filter Transformed and normalized count data to only contain Patch or Plaque and respective nonlesional controls

## Patch
Cases_Patch <- c("Patch_Pat3", "Patch_Pat5", "Patch_Pat8", "Patch_Pat10", "Patch_Pat13")
dummy_nonlesional <- paste("nonlesional", dummy, sep = "_")
dummy_Patch <- paste("Patch", dummy, sep = "_")
dummy <- c(dummy_nonlesional, dummy_Patch)
vsd_Patch <- select(
    vsd,
    all_of(dummy)
)

## Plaque
Cases_Plaque <- c("Plaque_Pat1", "Plaque_Pat4", "Plaque_Pat7", "Plaque_Pat9", "Plaque_Pat14")
dummy <- sub(".*_", "", Cases_Plaque)
dummy_nonlesional <- paste("nonlesional", dummy, sep = "_")
dummy_Plaque <- paste("Plaque", dummy, sep = "_")
dummy <- c(dummy_nonlesional, dummy_Plaque)
dummy <- dummy[-1]
vsd_Plaque <- select(
    vsd,
    all_of(dummy)
)

##################################
## Gene Set Enrichment Analysis ##
##################################

gsea_patch <-run_pathfindR(
    input = DESeq2_DEG_Patch,
    gene_sets = "costum",
    custom_genes = Reactome_PFR_interested,
    custom_descriptions = Reactome_PFR_interested_descriptions,
    visualize_enriched_terms = FALSE,
    plot_enrichment_chart = FALSE
)

gsea_plaque <-run_pathfindR(
    input = DESeq2_DEG_Plaque,
    gene_sets = "costum",
    custom_genes = Reactome_PFR_interested,
    custom_descriptions = Reactome_PFR_interested_descriptions,
    visualize_enriched_terms = FALSE,
    plot_enrichment_chart = FALSE
)

# Cluster the results
gsea_patch_clustered <- cluster_enriched_terms(
    enrichment_res = DESeq2_DEG_Patch,
    plot_clusters_graph = FALSE
)

gsea_plaque_clustered <- cluster_enriched_terms(
    enrichment_res = DESeq2_DEG_Plaque,
    plot_clusters_graph = FALSE
)

# Retain only the enriched pathways that are representitive for cluster
representative_Patch <- gsea_patch_clustered %>% dplyr::filter(Status == "Representative")
representative_Plaque <- gsea_plaque_clustered %>% dplyr::filter(Status == "Representative")


# Calculate Agglomerated and Z-Score Normalized Scores of Enriched Pathways for Each Sample
Score_Comparison_Patch <- score_terms(
    enrichment_table = representative_Patch,
    exp_mat = as.matrix(vsd_Patch),
    cases = Cases_Patch,
    plot_hmap = FALSE
)

Score_Comparison_Plaque <- score_terms(
    enrichment_table = representative_Plaque,
    exp_mat = as.matrix(vsd_Plaque),
    cases = Cases_Plaque,
    plot_hmap = FALSE
)


###########################
## Plot Score Comparison ##
###########################

# Define Colors for the Heatmap
col <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

##########
# Patch #
#########
score_matrix <- Score_Comparison_Patch
cases <- Cases_Patch
case_title = "Patch"
control_title = "nonlesional skin"

## Compute the Euclidean distance matrix
dist_matrix <- dist(score_matrix, method = "euclidean")

## Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "ward.D2")

## Reorder the rows based on the clustering results
score_matrix_reordered <- score_matrix[order.dendrogram(as.dendrogram(hc)), ]
score_matrix <- score_matrix_reordered 

## Create the heatmap using ggplot
tmp <- rowMeans(score_matrix[, cases, drop = FALSE])
score_matrix <- score_matrix[c(which(tmp >= 0), which(tmp < 0)), ]
var_names <- list()
var_names[["Term"]] <- factor(rownames(score_matrix), levels = rev(rownames(score_matrix)))
var_names[["Sample"]] <- factor(colnames(score_matrix), levels = colnames(score_matrix))
score_df <- expand.grid(var_names, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
scores <- as.vector(score_matrix)
scores <- data.frame(scores)
score_df <- cbind(score_df, scores)
score_df$Type <- ifelse(score_df$Sample %in% cases, case_title, control_title)
score_df$Type <- factor(score_df$Type, levels = c(case_title, control_title))


g <- ggplot2::ggplot(score_df, ggplot2::aes_(x = ~Sample, y = ~Term))
g <- g + ggplot2::geom_tile(ggplot2::aes_(fill = ~scores), color = "grey60")
g <- g + ggplot2::scale_fill_gradientn(
  colours = col,
  values = NULL,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "fill"
)
g <- g + ggplot2::theme(
  axis.title.x = ggplot2::element_blank(), 
  axis.title.y = ggplot2::element_blank(),
  axis.text.x = ggplot2::element_blank(),
  axis.ticks.x = ggplot2::element_blank(),
  legend.title = ggplot2::element_text(size = 10), 
  legend.text = ggplot2::element_text(size = 12)
)
g <- g + ggplot2::labs(fill = "Z-Score")
g <- g + ggplot2::facet_grid(~Type, scales = "free_x", space = "free")
g <- g + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, face = "bold"))

g <- g + theme(
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)


##########
# Plaque #
##########
score_matrix <- Score_Comparison_Plaque
cases <- Cases_Plaque
case_title = "Plaque"
control_title = "nonlesional skin"

## Compute the Euclidean distance matrix
dist_matrix <- dist(score_matrix, method = "euclidean")

## Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "ward.D2")

## Reorder the rows based on the clustering results
score_matrix_reordered <- score_matrix[order.dendrogram(as.dendrogram(hc)), ]
score_matrix <- score_matrix_reordered 

## Create the heatmap using ggplot
tmp <- rowMeans(score_matrix[, cases, drop = FALSE])
score_matrix <- score_matrix[c(which(tmp >= 0), which(tmp < 0)), ]
var_names <- list()
var_names[["Term"]] <- factor(rownames(score_matrix), levels = rev(rownames(score_matrix)))
var_names[["Sample"]] <- factor(colnames(score_matrix), levels = colnames(score_matrix))
score_df <- expand.grid(var_names, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
scores <- as.vector(score_matrix)
scores <- data.frame(scores)
score_df <- cbind(score_df, scores)
score_df$Type <- ifelse(score_df$Sample %in% cases, case_title, control_title)
score_df$Type <- factor(score_df$Type, levels = c(case_title, control_title))


g <- ggplot2::ggplot(score_df, ggplot2::aes_(x = ~Sample, y = ~Term))
g <- g + ggplot2::geom_tile(ggplot2::aes_(fill = ~scores), color = "grey60")
g <- g + ggplot2::scale_fill_gradientn(
  colours = col,
  values = NULL,
  space = "Lab",
  na.value = "grey50",
  guide = "colourbar",
  aesthetics = "fill"
)
g <- g + ggplot2::theme(
  axis.title.x = ggplot2::element_blank(), 
  axis.title.y = ggplot2::element_blank(),
  axis.text.x = ggplot2::element_blank(),
  axis.ticks.x = ggplot2::element_blank(),
  legend.title = ggplot2::element_text(size = 10), 
  legend.text = ggplot2::element_text(size = 12)
)
g <- g + ggplot2::labs(fill = "Z-Score")
g <- g + ggplot2::facet_grid(~Type, scales = "free_x", space = "free")
g <- g + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, face = "bold"))

g <- g + theme(
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)