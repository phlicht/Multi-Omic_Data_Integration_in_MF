library(ggplot2)
library(ggpubr)
library(EnhancedVolcano)
library(DESeq2)
library(dplyr)

###############
## Figure 1A ##
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
vsd <- vst(dds)

# Get PCA results
pcaData <- plotPCA(vsd, intgroup=c("Stage", "Patient"), returnData=TRUE)
pcaData$Stage <- factor(
    pcaData$Stage,
    levels = c(
        "nonlesional",
        "Patch",
        "Plaque"
    )
)
pcaData$Patient <- factor(
    pcaData$Patient,
    levels = c(
        "Pat1",
        "Pat3",
        "Pat4",
        "Pat5",
        "Pat7",
        "Pat8",
        "Pat9",
        "Pat10",
        "Pat13",
        "Pat14"
    )
)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot the PCA using ggplot2
## Define shapes and colors
Stage_shape <- c(
    "nonlesional" = 16,
    "Patch" = 17,
    "Plaque" = 15
)
Patient_Colour <- c(
    "Pat1" = "#F8766D",
    "Pat3" = "#D89000",
    "Pat4" = "#A3A500",
    "Pat5" = "#39B600",
    "Pat7" = "#00BF7D",
    "Pat8" = "#00BFC4",
    "Pat9" = "#00B0F6",
    "Pat10" = "#9590FF",
    "Pat13" = "#E76BF3",
    "Pat14" = "#FF62BC"
)


PCA_STAR_noSkin2 <- ggplot(
    pcaData,
    aes(PC1, PC2, color = Patient, shape = Stage)) +
    geom_point(size=4) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    scale_shape_manual(values = Stage_shape) +
    scale_color_manual(values = Patient_Colour) +
    theme_bw()

# Move legend to top
PCA_STAR_noSkin2 <- PCA_STAR_noSkin2 +
    guides(shape = guide_legend(direction = 'vertical')) +
    theme(
        legend.position = "top",
        legend.spacing.y = unit(-0.5, "line"),
        legend.title = element_blank(),
        legend.key.size = unit(1.2, units = "char"),
        legend.spacing = unit(1, "line"),
        legend.margin = margin(0,0,0,0),
        legend.box.spacing = unit(0, "pt")
    )

####################
## Figure 1B / 1C ##
####################

# Figure 1B
Patch_vulcano_plot <- EnhancedVolcano(
    dds_STAR_noSkin2_resPatchLFC_df_GeneSymbol,
    lab = rownames(dds_STAR_noSkin2_resPatchLFC_df_GeneSymbol),
    x = 'log2FoldChange',
    y = 'pvalue',
    pCutoffCol="padj",
    subtitle = NULL,
    caption = NULL,
    title = "Patch vs. nonlesional skin",
    legendPosition="right",
    legendLabSize = 10,
    legendIconSize = 2,
    titleLabSize = 12,
    labSize = 4 ,
    drawConnectors = T,
    axisLabSize = 14,
    captionLabSize = 12
)

# Figure 1C
Plaque_display_genes <- c(
    "KRT31",
    "KRT18",
    "CLDN10",
    "TP53AIP1",
    "NOP53",
    "TP53",
    "RASGRP1",
    "SLAMF7",
    "CD84",
    "MYCBP2",
    "ADAM17",
    "IFNG",
    "ICOS",
    "SERPINB3",
    "SERPINB4",
    "BIRC3",
    "PTPN22",
    "GPR15",
    "CD274"
)

Plaque_vulcano_plot <-EnhancedVolcano(
    dds_STAR_noSkin2_resPlaqueLFC_df_GeneSymbol,
    lab = rownames(dds_STAR_noSkin2_resPlaqueLFC_df_GeneSymbol),
    selectLab = Plaque_display_genes,
    x = 'log2FoldChange',
    y = 'pvalue',
    pCutoffCol = "padj",
    subtitle = NULL,
    caption = NULL,
    title = "Plaque vs. nonlesional skin",
    legendPosition = "right",
    legendLabSize = 10,
    legendIconSize = 2,
    titleLabSize = 12,
    labSize = 4,
    drawConnectors = T,
    axisLabSize = 14,
    captionLabSize = 12,
    pointSize = 1,
    max.overlaps = 23
)

## Arrange Figures 1B and 1C
Volcoano_arranged <- ggarrange(
    Patch_vulcano_plot,
    Plaque_vulcano_plot,
    common.legend = T,
    legend = "bottom",
    labels = c("B", "C")
)

##########################
## Arrange final Figure ##
##########################

figure_1 <- ggarrange(
    ggarrange(
        PCA_STAR_noSkin2,
        nrow = 1,
        ncol = 1,
        labels = "AUTO"
    ),
    Volcoano_arranged,
    nrow = 2
)

