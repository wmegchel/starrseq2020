####################################################################
#
# Boxplot of CpG-methylation in C1 and C2, WT and DNMT TKO cultured in SL
#
# Wout Megchelenbrink
# January 30, 2020
####################################################################
rm(list=ls())
source("include/style.R")

WT <- fread("DATA/CpG_DNA_Methylation_q10_rmdup_counts_SL_merged.tsv.gz")
TKO <- fread("DATA/CpG_methylation_SL_DNMT_TKO_GSM1891661.tsv.gz")

WT.gr <- makeGRangesFromDataFrame(WT, start.field = "position", end.field = "position")
TKO.gr <- makeGRangesFromDataFrame(TKO, start.field = "position", end.field = "position")


# Get the STARR-seq regions and concat random regions
starr   <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")[, .(chr, start=summit-250L, end=summit+250L, class)]
rnd     <- fread("DATA/RND_10_or_more_input_DNA_GC_matched.tsv.gz")
starr   <- rbind(starr, rnd)
starr[, ID:=.I]
starr.gr <- makeGRangesFromDataFrame(starr)

# Overlap STARR/RND peaks with CpG methylation
ovl.WT <- findOverlaps(starr.gr, WT.gr)
ovl.TKO <- findOverlaps(starr.gr, TKO.gr)

# Concat results
DT.WT <- cbind(starr[queryHits(ovl.WT), .(ID, chr, start, end, class)], WT[subjectHits(ovl.WT), .(nMeth, nTot, genotype="WT") ])
DT.TKO <- cbind(starr[queryHits(ovl.TKO), .(ID, chr, start, end, class)], TKO[subjectHits(ovl.TKO), .(nMeth, nTot, genotype="TKO") ])
DT <- rbind(DT.WT, DT.TKO)

# Take average CpG methylation per peak
DT <- DT[ , lapply(.SD, sum), by=.(ID, chr, start, end, class, genotype), .SDcols=c("nMeth", "nTot") ]
DT[, pct:=nMeth/nTot * 100]

DT[,genotype:=factor(genotype, levels = c("WT", "TKO"))]

# Make the barplot
ggplot(DT[class %in% c("C1", "C2")], aes(x=class, y=pct, color = genotype, fill = genotype)) + 
geom_boxplot(notch = F, outlier.shape =  NA,  cex = .25) +
stat_summary(geom = "crossbar", width=0.75, fatten=0, cex=.25, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
theme_wout_minimal() +
scale_fill_manual(values = colors.basic[c(4,6)]) +
scale_color_manual(values = colors.basic[c(4,6)]) +
xlab(NULL) +
ylab("CpG methylation (%)")

# Save PDF
ggsave("IMG/FigS3C_Barplot_CpG_methylation_C1_C2_WT_vs_KO.pdf", width = 2.5, height = 2)
