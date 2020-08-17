####################################################################################
#
# Barplot of the number of up- or downregulated genes 
# in P53KO compared to WT in 2iL and SL
#
# Wout Megchelenbrink
# January 22, 2020
#
####################################################################################

source("include/style.R")

# Read P53 WT/KO RNAseq
rnaseq <- fread("DATA/RNAseq_Trp53_WT_and_KO.tsv.gz")

# Get number and direction of DEGs
DT <- rbind(rnaseq[padj_KO_vs_WT_2iL < 0.05 & abs(log2FC_KO_vs_WT_2iL) >= log2(2.5), .(nGenes=.N, condition="2iL"), by=.(dir=ifelse(log2FC_KO_vs_WT_2iL < 0, "down", "up"))],
            rnaseq[padj_KO_vs_WT_SL < 0.05 & abs(log2FC_KO_vs_WT_SL) >= log2(2.5), .(nGenes=.N, condition="SL"), by=.(dir=ifelse(log2FC_KO_vs_WT_SL < 0, "down", "up"))])

DT[, dir.int:=ifelse(dir == "up", 1, -1)]

# Make the plot
ggplot(DT, aes(x=condition, y=nGenes*dir.int, fill = dir)) +
geom_bar(stat = "identity", color = "white") +
scale_fill_manual(values = colors.dark[c(2,3)], guide=F) +
scale_y_continuous(limits = c(-400,400)) +
theme_wout_minimal() +
xlab(NULL) +
ylab("Differential genes (#)")

# Save PDF
ggsave("IMG/Fig4E_P53_RNAseq_DEGs_up_down.pdf", width = 1.5, height = 2.5)  
