###################################################################################
#
# Scatterplot of STARR-seq signal at Zic3 ChIP-seq peaks in 2iL and SL
# DE peaks are purple
#
# Wout Megchelenbrink
# January, 20 2020
#
###################################################################################

rm(list=ls())
source("include/style.R")

FC.cutoff     <- 2.5 # STARR-seq peak with a fold change >= 2.5 are considered differential (if p < 0.05)
min.overlap   <- 250 # Consider STARR and ChIP peaks overlapping if they share >= 250 bp

# Read STARR-seq and annotate DE peaks
starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")[class %in% c("C1", "C2")]
starr[, isDE:=0]
starr[log2fc_SL_vs_2iL >= log2(FC.cutoff) & padj_SL_vs_2iL < .05, isDE:=1]

# Read ZIC3 ChIP-seq
zic3 <- fread("DATA/ZIC3_ChIP_and_STARRseq.tsv.gz")[ ,.(peak_id,chr,start,end,isPMT)]

# Overlap STARR-seq and ZIC3 ChIP-seq
starr.gr <- makeGRangesFromDataFrame(starr)
zic3.gr <- makeGRangesFromDataFrame(zic3)
ovl <- findOverlaps(starr.gr, zic3.gr, minoverlap = min.overlap)

# Assign the promoter/enhancer status based on the ZIC3 ChIP-seq peaks with 
# largest overlap
overlaps <- pintersect(starr.gr[queryHits(ovl)], zic3.gr[subjectHits(ovl)])
pctOverlap <- width(overlaps) / width(starr.gr[queryHits(ovl)]) * 100
DT <- cbind(starr[queryHits(ovl), ], zic3[subjectHits(ovl), .(peak_id, isPMT)], pctOverlap)
DT <- DT[DT[, .I[pctOverlap >= max(pctOverlap)], by=starr_peak_id]$V1,]
DT <- unique(DT, by="starr_peak_id")
DT[, isDE:=factor(isDE)]

# Split promoters and enhancers
DT[, isPMT:=ifelse(isPMT==1, "Promoter", "Enhancer")]
DT[, isPMT:=factor(isPMT, levels=c("Promoter", "Enhancer"))]

DTX <- DT[isDE == 1, .N, by=.(condition=ifelse(sign(log2fc_SL_vs_2iL)==1, "SL", "2iL"), isPMT)]
DTX <- rbind(DTX, data.table(condition=c("2iL", "2iL"),isPMT=c("Promoter","Enhancer"), N=c(0,0)))

ggplot(DTX, aes(x=condition, y=N)) +
geom_bar(stat = "identity", color = "white", fill = colors.basic[5]) + 
coord_flip() +
facet_wrap(~isPMT) +
xlab("") +
theme_wout_minimal()

ggsave("IMG/FigS2G_top_barplot_STARRseq_at_ZIC3_ChIP_peaks_2iL_vs_SL.pdf", width = 6, height = 1.5)


# Make the plot,
ggplot(DT, aes(x=log2(enrichment_2iL+1), y=log2(enrichment_SL+1), color = isDE)) +
geom_point(cex = .1) + 
coord_equal() +
scale_color_manual(values = c("#CCCCCC", colors.basic[5]), guide=F) +
theme_wout_minimal() +
scale_x_continuous(limits = c(0,6), breaks = seq(0,6,by=2)) +
scale_y_continuous(limits = c(0,6), breaks = seq(0,6,by=2)) +
geom_abline(slope = 1, linetype = "dashed", color = colors.basic[6]) +
xlab("STARR-seq 2iL (log2 enrichment)") +
ylab("STARR-seq SL (log2 enrichment)") +
facet_wrap(~isPMT)

# Save PDF
ggsave("IMG/FigS2G_bottom_scatterplot_STARRseq_at_ZIC3_ChIP_peaks_2iL_vs_SL.pdf", width = 6, height = 3)
