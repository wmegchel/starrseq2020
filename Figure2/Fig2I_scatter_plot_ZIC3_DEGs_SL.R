###################################################################################
#
# Scatterplot of DEGs in Zic3 KO vs WT
#
# Wout Megchelenbrink
# Feb, 14 2020
#
###################################################################################
source("include/style.R")
cutoff <- 250 # read count for expression

# Read data and convert to long format
DT <- fread("DATA/Zic3KO_RNAseq_2iL_SL.tsv.gz")[(WT_2iL1 >= cutoff | (KO_2iL1+KO_2iL2) / 2 >= cutoff) | (WT_SL1 >= cutoff | (KO_SL1+KO_SL2) / 2 >= cutoff)]
DTX <- melt.data.table(DT, id.vars = "gene_name", measure.vars = list(c("log2FC_2iL", "log2FC_SL"), c("pval_2iL", "pval_SL")), 
                       value.name = c("log2FC", "pval"), variable.name = "condition")
DTX[, condition:=ifelse(condition==1, "2iL","SL")]
DTX[, log2FC:=log2FC * -1]

# Annotate PrE markers
DTX[condition == "SL" & gene_name %in% c("Dab2", "Gata4", "Zic3", "Foxq1", "Sox17", "Cubn", "Lrp2", "Lama1", "Lamb1"), lab:=gene_name]
DTX[condition == "2iL" & gene_name %in% c("Zic3"), lab:=gene_name]
DTX[, col:=0L]
DTX[!is.na(lab), col:=1L]
DTX[, col:=factor(col)]

# Make the plot
ggplot(DTX[condition == "SL"], aes(x=log2FC, y=-log10(pval), label=lab, color=col)) +
geom_point_rast() + 
facet_wrap(~condition, nrow = 2) +
geom_text_repel(min.segment.length = 0) +
theme_wout_minimal() +
scale_x_continuous(limits = c(-10,10)) +
scale_y_continuous(limits = c(0,100)) +
scale_color_manual(values = c("#CCCCCC", colors.basic[2]))

# Save PDF
ggsave("IMG/Fig2I_scatterlot_ZIC3_DEGs.pdf", width = 3, height = 4)
