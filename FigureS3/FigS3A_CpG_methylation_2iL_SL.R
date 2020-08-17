####################################################################
#
# Wout Megchelenbrink
# CpG methylation in C1, C2, C3 and random regions
# Figure S3A
#
# Jan 30, 2020
####################################################################

rm(list=ls())
source("include/style.R")

# Read CpG data
WT.2iL <- fread("DATA/CpG_DNA_Methylation_q10_rmdup_counts_2iL_merged.tsv.gz")
WT.SL <- fread("DATA/CpG_DNA_Methylation_q10_rmdup_counts_SL_merged.tsv.gz")

# Make GRanges
WT.2iL.gr <- makeGRangesFromDataFrame(WT.2iL, start.field = "position", end.field = "position")
WT.SL.gr <- makeGRangesFromDataFrame(WT.SL, start.field = "position", end.field = "position")

# Get the STARR-seq regions and concat random regions
starr   <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")[, .(chr, start=summit-250L, end=summit+250L, class)]
rnd     <- fread("DATA/RND_10_or_more_input_DNA_GC_matched.tsv.gz")
starr   <- rbind(starr, rnd)
starr[, ID:=.I]
starr.gr <- makeGRangesFromDataFrame(starr)

ovl.2iL <- findOverlaps(starr.gr, WT.2iL.gr)
ovl.SL <- findOverlaps(starr.gr, WT.SL.gr)

DT.2iL <- cbind(starr[queryHits(ovl.2iL), .(ID, chr, start, end, class)], WT.2iL[subjectHits(ovl.2iL), .(nMeth, nTot) ])
DT.SL <- cbind(starr[queryHits(ovl.SL), .(ID, chr, start, end, class)], WT.SL[subjectHits(ovl.SL), .(nMeth, nTot) ])

DT.2iL[, culture:="2iL"]
DT.SL[, culture:="SL"]

DT <- rbind(DT.2iL, DT.SL)

# Sum CpG counts per pek and get fractions
DT <- DT[ , lapply(.SD, sum), by=.(ID, chr, start, end, class, culture), .SDcols=c("nMeth", "nTot") ]
DT[, pct:=nMeth/nTot * 100]

# Make the plot
ggplot(DT, aes(x=class, y=pct, fill = class, color = class)) + 
geom_boxplot(notch = F, outlier.shape =  NA,  cex = .25) +
stat_summary(geom = "crossbar", width=0.75, fatten=0, cex=.25, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
theme_wout_minimal() +
scale_color_manual(values = colors.basic[c(2,1,3,6)], guide=F) +
scale_fill_manual(values = colors.basic[c(2,1,3,6)], guide=F) +
xlab(NULL) +
ylab("CpG methylation (%)") +
facet_wrap(~culture)

# Save PDF
ggsave("IMG/Fig3A_CpG_methylation_C1_C2_C3.pdf", width = 2.5, height = 2)


# Stats
classes <- c("C1", "C2", "C3")
for(cond in c("2iL", "SL"))
{
  for(i in 1:length(classes))
  {
    for(j in 1:length(classes))
    {
      if(j > i)
      {
        cat(sprintf("Culture = %s , comparing %s vs %s :: p = %2.2E\n", cond, classes[i], classes[j],
                    wilcox.test(DT[culture == cond & class == classes[i], pct], DT[culture == cond & class==classes[j], pct])$p.value))
      }
    }
  }
}