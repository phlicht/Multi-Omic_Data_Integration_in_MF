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

###############
## Figure 1B ##
###############

# Patch volcano plot
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

# Plaque volcano plot
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

## Arrange both volcano plots
Volcoano_arranged <- ggarrange(
    Patch_vulcano_plot,
    Plaque_vulcano_plot,
    common.legend = T,
    legend = "bottom",
    labels = c("B")
)


###############
## Figure 1C ##
###############

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplotify)
library(ggpubr)
library(ggfortify)


## Preprocessing
################

# Select the pathways that may have been elicited by S. aureus (or other Staphylocci)
# according to pernitent literature
Saureus_pws_Plaque <- representative_Plaque %>% 
  dplyr::filter(
    Term_Description == "NOD1/2 Signaling Pathway" |
      Term_Description == "RHOA GTPase cycle" |
      Term_Description == "Interleukin-4 and Interleukin-13 signaling"
  ) %>%
  select(Term_Description, Up_regulated, Down_regulated)

# Extract genes of interest
up_regulated <- Saureus_pws_Plaque %>% 
  pull(Up_regulated) %>% 
  str_split(pattern = ", ")

names(up_regulated) <- Saureus_pws_Plaque %>% pull(Term_Description)


## Normalization and Transformation of Blood Samples
####################################################

# Since we want to assess differences between Plaque samples, we again use DESeq2's median of rations normalization
# and rlog transformation only on the raw gene count matrix of Plauqe samples 

# Filter metadata and row count matrix accordingly
metadata_Plaque <- metadata %>% 
  dplyr::filter(Stage == "Plaque") 

metadata_Plaque$Patient <- factor( metadata_Plaque$Patient, levels = unique(metadata_Plaque$Patient))

# Raw Counts
raw_counts <- raw_counts %>%
  t %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  dplyr::filter(Sample %in% metadata_Plaque$Sample) %>% 
  column_to_rownames("Sample") %>% 
  t

# Normalize and rlog transform
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = metadata,
  design= ~ Dominant_Species
)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
rld <- rlog(dds, blind = FALSE)


## Perform and Plot PCA
#######################
Patient_Colour <- c(
  "Pat1" = "#F8766D",
  "Pat4" = "#A3A500",
  "Pat7" = "#00BF7D",
  "Pat9" = "#00B0F6",
  "Pat14" = "#FF62BC"
)

# Perform PCA
pca_result <- mat %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL") %>%
  inner_join(ENSEMBL_to_HGNC) %>%
  dplyr::filter(Gene_symbol %in% unique(unlist(regulated))) %>% 
  column_to_rownames("Gene_symbol") %>% 
  select(!ENSEMBL) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  dplyr::filter(Sample %in% row.names(metadata_Plaque)) %>% 
  column_to_rownames("Sample") %>% 
  prcomp(center = TRUE, scale. = TRUE)


# Extract Data and process according to ggfortify::autoplot.prcomp
plot.data <- ggplot2::fortify(pca_result, data = metadata_Plaque)
plot.data$rownames <- rownames(plot.data)

ve <- object$sdev^2/sum(object$sdev^2)
PC <- paste0("PC", c(x, y))
x.column <- PC[1]
y.column <- PC[2]
loadings.column <- "rotation"
lam <- object$sdev[c(x, y)]
lam <- lam * sqrt(nrow(plot.data))

if (scale != 0) {
  lam <- lam^scale
  plot.data[, c(x.column, y.column)] <- t(t(plot.data[, c(x.column, y.column)])/lam)
}
plot.columns <- unique(c(x.column, y.column, colnames(plot.data)))
plot.data <- plot.data[, plot.columns]
if (!is.null(loadings.column)) {
  loadings.data <- as.data.frame(object[[loadings.column]][, 
  ])
  loadings.data$rownames <- rownames(loadings.data)
  loadings.columns <- unique(c(x.column, y.column, colnames(loadings.data)))
  loadings.data <- loadings.data[, loadings.columns]
} else {
  loadings.data <- NULL
}
if (is.null(ve) | !variance_percentage) {
  labs <- PC
} else {
  ve <- ve[c(x, y)]
  labs <- paste0(PC, " (", round(ve * 100, 2), "%)")
}
xlab <- labs[1]
ylab <- labs[2]

ggplot(plot.data, aes(x = PC1, y = PC2, color = Patient, shape = Dominant_Species)) +
  geom_point(size = 4, alpha = 0.8, position = position_jitter(width = 0.013, height = 0.013)) +
  labs(
    x = xlab,
    y = ylab,
    shape = "Dominant Species",
    caption = paste0(length(unique(unlist(regulated))), " genes:\nNOD1/2 Signalling, RHOA GTPase cycle, IL-13 Signalling")
  ) +
  ggtitle(
    "PCA of genes in Pathways that may have been elicited \nby S. aureus colonization",
    subtitle = "Plaque Samples"
  ) +
  scale_color_manual(values = Patient_Colour) +
  theme_bw() 