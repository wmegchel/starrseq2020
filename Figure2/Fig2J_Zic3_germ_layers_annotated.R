##################################################################
#
# Heaatmap that shows the expression of some typical developmental 
# marker genes in  WT and ZIC3 KO ESCs
#
# Wout Megchelenbrink
#
# Feb, 14 2020
#
##################################################################
library(circlize)
rm(list=ls())
source("include/style.R")

# Read the data
cutoff <- 250
DT <- fread("DATA/Zic3KO_RNAseq_2iL_SL.tsv.gz")[(WT_2iL1 >= cutoff | (KO_2iL1+KO_2iL2) / 2 >= cutoff) | (WT_SL1 >= cutoff | (KO_SL1+KO_SL2) / 2 >= cutoff)]
x <- DT[order(padj_SL)][log2FC_SL > 0][1:200]

# Define cluster markers / developmental genes
pluri     <- c("Pou5f1", "Sox2", "Nanog", "Klf2", "Klf4", "Klf5", "Esrrb", "Zfp42", "Prdm14", "Tcfcp2l1", "Zic3")
epiblast  <- c("Lefty1", "Lefty2")
endo      <- c("Gata4", "Gata6", "Sox7", "Sox17", "Hnf4a", "Foxa2", "Dab2", "Foxq1", "Cubn", "Lrp2", "Lama1", "Lamb1")

# Assign genes to clusers
DT[gene_name %in% pluri, type:="Pluri"]
DT[gene_name %in% epiblast, type:="Epi"]
DT[gene_name %in% endo, type:="Endo"]

DT <- DT[!is.na(type)]
DT[type == "Pluri"]


# Log2 transform
mat <- log2(data.matrix(DT[, .(WT_2iL1, KO_2iL1, KO_2iL2, WT_SL1, KO_SL1, KO_SL2)]) + 1)
rownames(mat) <- DT$gene_name

# Clusters
clusters <- c("Pluri", "Endo", "Epi")

# Get the ones that are specifc
partition <- DT[, .(type)]
partition[, type:=factor(type, levels = clusters)]

# Cluster annotation
hm.clust <- Heatmap(data.matrix(partition$type), 
                    col = structure(colors.basic, names = clusters), 
                    name = "Cluster", show_row_names = F, width = unit(4, "mm"), show_heatmap_legend = F)


# Make the heatmap
col.max <- 15 # max color intensity = 2^15 normalized read counts
hm.list <- Heatmap(mat, clustering_distance_rows = "euclidean" , cluster_columns = F, show_row_dend = F, show_row_names = T, row_names_side = "left",
                   col = colorRamp2(c(0, 7.5, col.max), c("#1600aa", "white", "#c60000")), gap = unit(3, "mm"),
                   heatmap_legend_param = list(title = "RNAseq Zic3-/- vs (log2 read count)", color_bar = "continuous", legend_direction = "horizontal",
                                               legend_width = unit(50, "mm"), grid_border = "#555555", at = c(0, 7.5, col.max),
                                               labels = c(0, 7.5, sprintf("> %2.1f", col.max)))) + hm.clust
 
# Save PDF
pdf("IMG/Fig2J_Heatmap_RNAseq_Zic3_KO_germ_layers.pdf", width = 4, height = 7)
ComplexHeatmap::draw(hm.list,   heatmap_legend_side = "bottom", padding = unit(c(5,5,5,5), "mm") ,  split = partition)
dev.off()