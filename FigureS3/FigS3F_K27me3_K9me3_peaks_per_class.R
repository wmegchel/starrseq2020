####################################################################
#
# Wout Megchelenbrink
# H3K9me3, H3K27me3 peaks overlapping C1, C2, C3 and random peaas
#
# Jan 30, 2020
####################################################################
rm(list=ls())
source("include/style.R")

# Get STARR and random classes 
starr<- fread("DATA/STARRseq_OSN_RND_250bp.tsv.gz")
starr.gr  <- makeGRangesFromDataFrame(starr)

# Get K27me3, K9me3 peaks
DT <- fread("DATA/H3K9me3_H3K27me3_broad_peaks_P12_B5.tsv.gz") 
DT.gr <- makeGRangesFromDataFrame(DT)

# Overlap and concat
ovl <- findOverlaps(starr.gr, DT.gr)
DTX <- cbind(starr[queryHits(ovl), ], DT[subjectHits(ovl), ])

DTX <- unique(DTX, by=c("starr_peak_id","mark","condition"))
DTY1 <- merge(DTX[, .(N=length(unique(starr_peak_id ))), by=.(mark, condition, Class)], starr[, .(Ntot=.N), by=Class], by="Class")
DTY1[, pct:=N/Ntot * 100]

DTY2 <- data.table(Class=c("C1","C2","C3","RND"), mark=rep("H3K27me3",4), condition=rep("2iL",4), N=rep(0,4), Ntot=c(7072,18544,1950,10000), pct=rep(0,4))
DTY <- rbind(DTY1, DTY2)

# Make the plot
ggplot(DTY, aes(x=Class, y=pct, fill=condition)) +
geom_bar(stat = "identity", position = "dodge", color = "white") +
facet_wrap(~mark) +
scale_y_continuous(limits = c(0,20), breaks=c(0,10,20)) +
scale_fill_manual(values = colors.basic[1:2]) +
coord_flip() +
theme_wout_minimal()

# Save PDF
ggsave("IMG/FigS3F_Barplot_K27me3_K9me3_per_class.pdf", width=4, height = 3)
