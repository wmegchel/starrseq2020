####################################################################
#
# author:   Wout Megchelenbrink
# sumary:   Pie chart of the genomic location of STARR-seq peaks
#
# Januari 20, 2020
####################################################################

rm(list=ls())
source("include/style.R")

# Read STARR-seq data
starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")
starr[, feature:=factor(feature, levels=c("Promoter", "Exon", "Intron", "Intergenic"))]

starr[ , in2iL:=F]
starr[enrichment_2iL >= 3 & padj_2iL1 < 0.05 & padj_2iL2 < 0.05, in2iL:=T]

starr[ , inSL:=F]
starr[enrichment_SL >= 3 & padj_SL1 < 0.05 & padj_SL2 < 0.05, inSL:=T]

# Count fraction of peaks per feature
starr.2iL <- starr[in2iL==T, .N, by=feature]
starr.2iL[ ,frac:=N/sum(N)]
starr.2iL[, sum(N)]

# Make the plot
ggplot(starr.2iL, aes(x="", y=N, fill=feature)) + 
geom_bar(width = 1, stat = "identity", color="white") +
coord_polar("y", start=0) +
scale_fill_manual(values = colors.dark[4:1]) +
theme(axis.text.x=element_blank()) +
geom_text(aes(y = N,label = percent(frac)), size=5) +
theme_void()

# Save PDF
ggsave("IMG/Fig1B_piechart_STARRseq_peaks_per_genomic_feature_2iL.pdf", width = 3.5, height = 3.5)


# Count fraction of peaks per feature
starr.SL <- starr[inSL==T, .N, by=feature]
starr.SL[ ,frac:=N/sum(N)]
starr.SL[, sum(N)]

# Make the plot
ggplot(starr.SL, aes(x="", y=N, fill=feature)) + 
geom_bar(width = 1, stat = "identity", color="white") +
coord_polar("y", start=0) +
scale_fill_manual(values = colors.dark[4:1]) +
theme(axis.text.x=element_blank()) +
geom_text(aes(y = N,label = percent(frac)), size=5) +
theme_void()

# Save PDF
ggsave("IMG/Fig1B_piechart_STARRseq_peaks_per_genomic_feature_SL.pdf", width = 3.5, height = 3.5)
