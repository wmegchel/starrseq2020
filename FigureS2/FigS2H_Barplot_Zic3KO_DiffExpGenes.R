#######################################################################################
#
# Barplots of DEGs in ZIC3 KO vs WT (up- or down)
#
# Wout Megchelenbrink
# January 23, 2020
#
#######################################################################################
source("include/style.R")

# Read Zic3 WT and KO RNA-seq
DT <- fread("DATA/Zic3KO_RNAseq_2iL_SL.tsv.gz")

# Convert to long
DTX <- melt.data.table(DT, 
                       id.vars = "gene_name", 
                       measure.vars = list(log2FC=c("log2FC_2iL", "log2FC_SL"), padj=c("padj_2iL", "padj_SL")),
                       variable.name = c("condition"))
DTX[, log2FC:=log2FC * -1]
# Take DE and annotate direction and condition
DTX[, condition:=ifelse(condition==1, "2iL", "SL") ]
DTX <- DTX[abs(log2FC) >= log2(2.5) & padj < .05, .N, by=.(sgn=sign(log2FC), condition)]
DTX[, dir:=factor(ifelse(sgn == 1, "upregulated", "downregulated"), levels = c("upregulated", "downregulated"))]

# Make the plot
ggplot(DTX, aes(x=condition, y=N*sgn, fill=dir)) +
geom_bar(stat="identity", color = "white") +
scale_fill_manual(values = colors.dark[3:2]) +
theme_wout_minimal() +
theme(legend.title = element_blank()) +
xlab(NULL) +
ylab("DE Genes (#)")

# Save PDF
ggsave("IMG/Barplot_ZIC3_DEGs.pdf", width = 4, height = 6)
