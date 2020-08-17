###########################################################################################
# scRNA-seq UMAP and features per cluster. Somehow clusters differently on MAC and Ubuntu,
# though Seurat is version 3.0.0 on both. UMAP version different?
# 
# Wout Megchelenbrink
#
# Last update: January 28, 2020
###########################################################################################
rm(list=ls())

library(Seurat)

source("include/style.R")
# Load RDS file
Zic3KO <- readRDS("DATA/scRNAseq_Zic3_WT_and_KO.Rds")

# Subset
Zic3KO <- subset(x = Zic3KO, subset = nCount_RNAseq >= 10000 &
                                      nFeature_RNAseq >= 2000 & 
                                      percent.mito < 25)
Zic3KO <- NormalizeData(object = Zic3KO, normalization.method = "LogNormalize", scale.factor = 10000)
Zic3KO <- ScaleData(object = Zic3KO)

Zic3KO <- FindVariableFeatures(Zic3KO)
Zic3KO <- RunPCA(object = Zic3KO)

Zic3KO <- RunUMAP(Zic3KO, dims =1:12)
UMAPPlot(Zic3KO, group.by = "plate", shape = "medium")
ggsave(filename = "IMG/Fig2I_top_UMAP_Zic3KO_12dims.svg", width = 7, height = 4)

Zic3KO <- FindNeighbors(Zic3KO, dims = 1:12)
Zic3KO <- FindClusters(Zic3KO, resolution = 0.1)
UMAPPlot(object = Zic3KO,  shape = "medium") +
scale_color_manual(values = colors.basic)
ggsave(filename = "IMG/Fig2I_bottom_UMAP_Zic3KO_12dims_clustered_0.1_resolution_mod_0.95.svg", width = 7, height = 4)

DEGs <- as.data.table(FindMarkers(Zic3KO, ident.1 = "1", ident.2 = "2"), keep.rownames ="gene")
DEGs[gene %in% c("Nanog", "Tcfcp2l1", "Klf2", "Prdm14", "Esrrb")]

### Heatmap for selected cluster markers #####
genes <- c("Nanog", "Tcfcp2l1", "Tbx3", "Prdm14", "Tet2", # 2iL WT
         "Zic3",  "Lefty1",  "Skil",  "Sox2", "Emb", "Trh", "Gbx2", "Zfp462", # SL WT
         "Trim28", "Dnmt3a", "Dnmt3b", "Dnmt3l", "Bmp4", "Utf1",  # SL KO
          "Sox7", "Sox17", "Gata4", "Gata6") # PrE

split <- c(rep(1, 5),  rep(2,8), rep(3, 6), rep(4,4))

# Grep these guys from the expression matrix; annotate genotype, cluster, medium; do hierarch clustering
mat <- Seurat::GetAssayData(Zic3KO, assay = "RNAseq", slot = "scale.data")
mat <- mat[genes, ]

cols.cluster <- colors.basic[c(1:5)]
names(cols.cluster) <- 0:4

cols.genotype <- c(colors.basic[4], "#666666")
names(cols.genotype) <- c("WT", "Zic3KO")

cols.medium <- c("#222222", "#BBBBBB")
names(cols.medium) <- c("ser", "2i")


Zic3KO$seurat_clusters <- factor(Zic3KO$seurat_clusters, levels = c(0,2,1,3))

# Order clusters
ID <- c(which(Zic3KO$seurat_clusters == 0),
        which(Zic3KO$seurat_clusters == 2),
        which(Zic3KO$seurat_clusters == 1),
        which(Zic3KO$seurat_clusters == 3))
length(which(Zic3KO$seurat_clusters == 1))
norm.data = Zic3KO@assays$RNAseq@counts

# Annotate clusters
ha <- HeatmapAnnotation(cluster = factor(Zic3KO$seurat_clusters[ID]),
                        medium = Zic3KO$medium[ID],
                        genotype = Zic3KO$genotype[ID],
                        col = list(cluster = cols.cluster,
                                   genotype = cols.genotype,
                                   medium = cols.medium),
                                   gp = gpar(border=F, col=NA, cex=0, lwd=0))

# Make heatmap and save as PDF
pdf("IMG/Fig2K_Heatmap_clusters_0_4.pdf", width = 5, height = 4)
Heatmap(mat[, ID], name = "Heatmap", cluster_rows = F, cluster_columns = F, col = colorRamp2(seq(-2,2,length.out = 10), PurpleAndYellow(10)),
        show_column_names = F, show_row_dend = F, show_column_dend = F, top_annotation = ha, split=split,  rect_gp = gpar(border=F, col=NA, cex=0, lwd=0), 
        use_raster = T, raster_quality = 2, raster_device = "CairoPNG")
dev.off()
