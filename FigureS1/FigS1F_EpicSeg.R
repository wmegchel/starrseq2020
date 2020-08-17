#######################################################################################
#
#
# EpiCSeg cluster map
# Wout Megchelenbrink
#
# January 23, 2020
#######################################################################################
rm(list=ls())
library(circlize)
source("include/style.R")

STARR <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")
epi.2iL <- STARR[, .N, by=epicseg_2iL][order(epicseg_2iL), N]
epi.SL <- STARR[, .N, by=epicseg_SL][order(epicseg_SL), N]

# Read the data and convert to matrix
epi <- fread("DATA/EpiCSeg_5means.tsv.gz", header = T)
mat <- t(data.matrix(epi[, 2:ncol(epi)]))
colnames(mat) <- epi$mark

# Convert to Z-cores
mat <- scale(mat, center =T)
rownames(mat) <- c("Bivalent", "Enhancer", "Promoter", "Heterochromatin", "Empty")

# Truncate negative values
mat[mat < 0] <- 0

row_ha = rowAnnotation(bar.2iL = anno_barplot(epi.2iL), bar.SL=anno_barplot(epi.SL))


# Make heatmap
hm <- Heatmap(mat, name ="Z+ score", col = colorRamp2(c(0, 2), c("white", colors.basic[2])), cluster_rows = F,  #rect_gp = chm.border.style.dark, 
              row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), show_column_dend =T,
              heatmap_legend_param = list(at=c(0,2), legend_height = unit(3, "cm")),
              right_annotation = row_ha) 

# Save PDF
pdf("IMG/FigS1F_Episeg_5states.pdf", width = 6.5, height = 2.25)
draw(hm)
dev.off()
