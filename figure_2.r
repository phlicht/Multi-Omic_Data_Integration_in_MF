library(ggplot2)
library(ggpubr)
library(pheatmap)
library(MOFA2)
library(DESeq2)
library(dplyr)

# Read in genes of interest to plot
genes_of_interest <- readr::read_tsv("genes_of_interest.tsv")

# Define Annotation Colours
annotation_colour <- list(
    Stage = c(
        nonlesional = "#0000CD",
        Patch = "#8F9C1B",
        Plaque = "#2D5A27"
    )
)

###############
## Figure 2A ##
###############

# Load Data into DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = metadata,
    design = ~Patient+Stage
)

# Normalize using DESeq2
dds <- DESeq2(dds)

# Transform DESeq2 normalized data using variance stabilizing transformation
vsd <- rlog(dds)
vsd <- assay(vsd)

# Filter input data.frame to genes of interest
df <- subset(
    vsd,
    rownames(vsd) %in% genes_of_interest_presentInMOFA$ENSEMBL
) 

# Plot
heatmap_not_denoised <- pheatmap(
    df,
    show_rownames = FALSE,
    show_colnames = FALSE, 
    treeheight_row = 0,
    treeheight_col = 20,
    annotation_col  = metadata[1], # Filtering metadata to only annotate Stage
    annotation_colors = annotation_colour,
    main = "Transcriptome"
)



###############
## Figure 2B ##
###############
heatmap_denoised <- plot_data_heatmap(
    MOFA2_model.trained,
    view = "Transcriptome",
    factor = 4,
    features = genes_of_interest$ENSEMBL,
    show_rownames = FALSE,
    show_colnames = FALSE,
    denoise = TRUE,
    treeheight_row = 0,
    treeheight_col = 20,
    annotation_samples = metadata[1], # Filtering metadata to only annotate Stage
    annotation_colors = annotation_colour,
    main = "Transcriptome (+ Microbiome Integration)"
)

##########################
## Arrange final Figure ##
##########################

figure_2 <- ggarrange(
    ggplotify::as.ggplot(heatmap_not_denoised),
    ggplotify::as.ggplot(heatmap_denoised),
    legend = "bottom"
)