########################################################################################################
#
# Nucleosome free regions at the 500bp center of P53-C1, P53-C2, OSN and random regions
#
# Wout Megchelenbrink
# Jan 30, 2020
########################################################################################################
rm(list=ls())
source("include/style.R")

dist.max <- 1e9

p53     <- fread("DATA/P53_ChIPseq_and_STARRseq.tsv.gz")[class_2iL == class_SL]
p53.gr  <- makeGRangesFromDataFrame(p53, seqnames.field = "chr", start.field = "motif_start", end.field = "motif_end")

target.genes  <- fread("DATA/RNAseq_Trp53_WT_and_KO.tsv.gz")[(log2FC_KO_vs_WT_2iL <= -log2(2.5) & padj_KO_vs_WT_2iL < .05) | (log2FC_KO_vs_WT_SL <= -log2(2.5) & padj_KO_vs_WT_SL < .05), unique(gene_name) ]
pmt           <- fread("DATA/Promoters.tsv.gz")[ ,.(transcript_id, gene_name, chr, start=tss, end=tss)][gene_name %in% target.genes]
pmt.gr        <- makeGRangesFromDataFrame(pmt)

# Find overlaps
ovl <- findOverlaps(pmt.gr, p53.gr, maxgap = dist.max)
DT  <- cbind(pmt[queryHits(ovl), .(gene_name, start,end) ], p53[subjectHits(ovl), .(peak_id, motif_start, motif_end, class=class_2iL)])
DT[, dist:=abs((motif_start+motif_end)/2 - (start+end)/2) ]

# Take the closest distance to the gene
DT <- DT[DT[, .I[dist <= min(dist)], by=.(gene_name, class)]$V1, ]
DT <- unique(DT, by=c("gene_name", "class"))            
DT[, dist:=log10(dist+1)]

# Make the plot
ggplot(DT, aes(x=class, y=dist,  color = class)) +
geom_violin() +
theme_wout_minimal() +
scale_color_manual(values = colors.basic[c(2,1)], guide=F) +
scale_y_continuous(limits = c(0,7), breaks = c(0,2,4,6)) +
ylab("Distance (log10)") +
theme_wout_minimal()

# Save PDF
ggsave("IMG/FigS4E_Distance_TSS_downreg_gene_to_P53_peak.pdf", width = 2, height = 3.5)

# Stats
wilcox.test(DT[class == "C1", dist], DT[class == "C2", dist])$p.value
