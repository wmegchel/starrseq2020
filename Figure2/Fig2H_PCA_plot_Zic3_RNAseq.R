###################################################################################
#
# PCA plot of RNA-seq profiles in Zic3 WT and KO ESCs. Top 500 variable genes
#
# Wout Megchelenbrink
# Feb, 14 2020
#
###################################################################################
source("include/style.R")

vst <- readRDS("DATA/Zic3_RNAseq_vst_transformed.Rds")
DT <- DESeq2::plotPCA(vst, intgroup = c("genotype", "condition"), ntop = 500, returnData = T)

source("include/style.R")
ggplot(DT, aes(x=PC1, y=PC2, color=condition, shape=genotype)) +
geom_point(cex = 2) +
scale_shape_manual(values = c(1,6)) +
scale_color_manual(values = colors.basic[1:4]) +
theme_wout_minimal() +
xlab("PC1: 83% variance") +
ylab("PC2: 13% variance") +
theme(legend.title = element_blank())

ggsave("IMG/Fig2H_PCA_Zic3KO.svg", width = 3.5, height = 2)  