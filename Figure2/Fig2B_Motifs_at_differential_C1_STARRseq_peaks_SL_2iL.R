###################################################################################
#
# Scatterplot of STARR-seq signal in 2iL and SL, with DE peaks annotated
#
# Wout Megchelenbrink
# January, 20 2020
#
###################################################################################

rm(list=ls())

source("include/style.R")

# Read the Homer2 motif data
DT <- fread("DATA/Homer_Motifs_C1_DE_STARRseq_peaks.tsv.gz")[!(abs(log2_enrichment) > 1.5 & -log10(pval) < 5)]
DT[, score:=-log10(padj)]
DT[class == "C1_SL", log2_enrichment:=log2_enrichment*-1]

# Annotate some outliers
DT[score > 40, TF:=motif]
DT[score > 20 & abs(log2_enrichment) >= 1, TF:=motif]

# Make the plot
# Unknown-ESCs element and ZIC are the same as ZIC3 (almost the same motif)
ggplot(DT, aes(x=log2_enrichment, y=score, label=TF)) +
geom_point(color = "#CCCCCC") + 
geom_text_repel() +
xlab("TF motif enrichment (log2 FC)") +
ylab("Enrichment (-log10([p-value])") +
theme_wout_minimal() 

# Save PDF
ggsave("C2/Motif_Analysis/img/C2_Motifs_SL_2iL.svg", width = 4, height = 3.5)