########################################################################
# 
# Barplots showing that H3K27ac is transiently gained at P53-C2
# sites during iPSC repgramming
#
# Wout Megchelenbrink
# Janary 22, 2020
#
########################################################################

rm(list=ls())
source("include/style.R")

k27ac.peaks <- fread("DATA/H3K27ac_iPSC_repgramming_Zviran_Chen_P15_B7.tsv.gz")
k27ac.peaks.gr <- makeGRangesFromDataFrame(k27ac.peaks)

P53.C2 <- fread("DATA/P53_OSN_RND.saf.gz")[, .(peak_id, chr, start=start_250, end=end_250, class)][class == "C2"]
P53.C2[, Ntot:=.N, by=class]
P53.C2.gr <- makeGRangesFromDataFrame(P53.C2)

ovl <- findOverlaps(P53.C2.gr, k27ac.peaks.gr)
DT <- cbind(P53.C2[queryHits(ovl), .(peak_id, class, Ntot)], k27ac.peaks[subjectHits(ovl), .(sample, dataset)])

DTX <- unique(DT[, .(frac.k27ac=length(unique(peak_id))/Ntot), by=.(class, sample, dataset)], by=c("class", "sample", "dataset"))


######## Plot the data for the Zviran et al (Mbd3 F/-) iPSC repgramming data set ########
zviran <- DTX[dataset == "zviran"]
zviran[, ID:=str_sub(sample, end=2)] 
setorder(zviran, ID)
zviran[, sample:=factor(sample, levels=unique(sample))]

ggplot(zviran, aes(x=sample, y=frac.k27ac)) +
geom_bar(stat="identity", position = "dodge", color = "white", fill = colors.basic[1]) +
theme_wout_minimal() +
theme(axis.text.x = element_text(angle = 90)) +
xlab(NULL) +
ylab("P53-C2 loci overlapping H3K27ac peak") +
scale_y_continuous(labels=percent, limits = c(0,.6), breaks = c(0,.2,.4,.6))

# Save PDF
ggsave("IMG/Fig4H_Barplot_Zviran_iPSC_H3K27ac_at_P53_C2.pdf", width = 10, height = 5)

######## Plot the data for the Chen et al iPSC repgramming data set ########
chen <- DTX[dataset == "chen"]
chen[, ID:=str_sub(sample, end=2)]
setorder(chen, ID)
chen[, sample:=factor(sample, levels=unique(sample))]

ggplot(chen, aes(x=sample, y=frac.k27ac)) +
geom_bar(stat="identity", position = "dodge", color = "white", fill = colors.basic[1]) +
scale_fill_manual(values = colors.basic[1]) +
xlab(NULL) +
ylab("P53-C2 loci overlapping H3K27ac peak") +
scale_y_continuous(labels=percent, limits = c(0,.6), breaks = c(0,.2,.4,.6)) +
theme_wout_minimal() + 
theme(axis.text.x = element_text(angle = 90)) 

# Save PDF  
ggsave("IMG/Fig4H_Barplot_Chen_iPSC_H3K27ac_at_P53_C2.pdf", width = 10, height = 5)
