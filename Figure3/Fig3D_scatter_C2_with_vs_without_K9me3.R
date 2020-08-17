####################################################################
#
# Motif enrichment for NRF1 motifs with/without H3K9me3
# Wout Megchelenbrink
#
# January 22, 2020
####################################################################

rm(list=ls())

source("include/style.R")

DT <- fread("DATA/Homer_motif_enrichment_C2_with_or_without_H3K9me3.tsv.gz")[log2_enrichment >= log2(1.5) & pval < 1e-10 & pct.target >= 3][order(pval)]
DT[grep("yy1", motif, ignore.case = T), TF:="YY1"]

# Remove motifs that are from the same 
DT <- DT[str_length(TF) > 0]
DT[, log10p:=-log10(pval)]
DT[, mp:=ifelse(class == "with_K9me3", 1, -1)]

# Make the plot
ggplot(DT, aes(x=2^log2_enrichment, y=log10p*mp, label=TF)) +
geom_point() +
theme_wout_minimal() +
geom_text_repel() +
scale_x_continuous(limits = c(0,4), breaks = c(0,2,4)) +
xlab("Fold enrichment") +
ylab("Significance (-log10[p-value])") +
geom_hline(yintercept = 0, linetype = "dashed", color = colors.basic[6])
  
# Save PDF
ggsave("IMG/Fig3D_C2_motif_enrichment_with_vs_without_K9me3.pdf", width = 3, height = 3)
