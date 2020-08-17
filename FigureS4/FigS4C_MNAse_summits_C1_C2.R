########################################################################
# 
# MNAse-seq summit occupancy and fuziness computed by Danpos2
# Wout Megchelenbrink
# 
# January 30, 2020
#
########################################################################
rm(list=ls())
max.dist.peakcenter.to.summit <- 74

source("include/style.R")

peaks <- fread("DATA/P53_OSN_RND.saf.gz")[class != "SWITCH"]
peaks.gr <- makeGRangesFromDataFrame(peaks, seqnames.field = "chr", start.field = "start_250", end = "end_250")
peaks[, .N, by=class]

MNAse <- fread("DATA/Danpos2_MNAse_SL.tsv.gz")
MNAse.gr <- makeGRangesFromDataFrame(MNAse)
ovl <- findOverlaps(peaks.gr, MNAse.gr)

DT <- cbind(peaks[queryHits(ovl), .(peak_id, chr, start_250, end_250, class) ], MNAse[subjectHits(ovl), .(smt_pos, smt_value, fuzziness_score)])
DT[, dist:=abs(round((start_250 + end_250)/2) - smt_pos)]
DT <- DT[dist <= max.dist.peakcenter.to.summit]

# Take the summit of the nucleosome closest to the peak center
DT <- DT[DT[, .I[dist <= min(dist)], by=peak_id]$V1,]

# In the rare case that two summits are eqal distant from the peak center, take the highest nucpos summit
DT <- DT[DT[, .I[smt_value >= max(smt_value)], by=peak_id]$V1,]
DT <- unique(DT, by="peak_id")
DT[, class:=factor(class, levels=c("C1", "C2", "OSN", "RND"))]

# Make the plot
ggplot(DT, aes(x=class, y=log2(smt_value), color=class, fill = class)) +
geom_boxplot(position = position_dodge(), outlier.shape = NA, notch=T) +
stat_summary(geom = "crossbar",position = position_dodge(), fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
scale_fill_manual(values = colors.dark[c(1,2,4,6)], guide=F) +
scale_color_manual(values = colors.dark[c(1,2,4,6)], guide=F) +
theme_wout_minimal() +
coord_cartesian(ylim = c(5.9,10.5)) +
xlab("") +
ylab("Nucleosome summit occupancy (log2 reads)")

# Save PDF
ggsave("IMG/FigS4C_Danpos_NucOccupancy_summit_74bp.pdf", width=3, height = 2.5)


# Get p-values
classes <- c("C1", "C2", "OSN", "RND") 
for(i in 1:length(classes) )
{
  for(j in 1:length(classes))
  {
    if(i < j)
      cat(sprintf("%s vs %s => p = %2.2E\n", classes[i], classes[j], wilcox.test(DT[class == classes[i], log2(smt_value+1)], DT[class == classes[j], log2(smt_value+1)])$p.value))    
  }
}

# Summary stats
DT[, .(avg=mean(smt_value), fuz=mean(fuzziness_score), N=length(unique(peak_id))), by=class]
