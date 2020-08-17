####################################################################
#
# Overlap of P53-C1, P53-C2, OSN and random peaks with ATAC-seq (2iL/SL)
#
# Wout Megchelenbrink
# Jan 30, 2020
####################################################################

rm(list=ls())
source("include/style.R")

# Read C1, C2, OSN and RND peaks
peaks <- fread("DATA/P53_OSN_RND.saf.gz")[class != "SWITCH"]
peaks[, center:=as.integer((start_250 + end_250)/2)]          

# Read ATAC-seq peaks
atac <- fread("DATA/ATACseq_IDR_0.03_2iL_and_SL.tsv.gz")
atac.gr <- makeGRangesFromDataFrame(atac)
peaks.gr <- makeGRangesFromDataFrame(peaks, start.field = "center", end.field = "center")

# Overlap
ovl <- findOverlaps(peaks.gr, atac.gr)

# Concat
DT <- cbind(peaks[queryHits(ovl), .(peak_id, chr, start, end, class) ], atac[subjectHits(ovl), .(condition)])

# Get statistics
stat <- merge(DT[, .(Novl=length(unique(peak_id))), by=.(class, condition)], peaks[, .N, by=.(class)], by="class")
stat[, pct:=Novl/N ]
stat[, class:=factor(class, levels=c("RND","OSN", "C2", "C1"))]

# Make barplot
ggplot(stat, aes(x=class, y=pct, fill=class)) +
geom_bar(stat="identity", position = "dodge") +
scale_fill_manual(values = colors.basic[rev(c(2,1,4,6))], guide=F) +
coord_flip() +
facet_wrap(~condition, nrow = 2) +
theme_wout_minimal() +
xlab("") +
scale_y_continuous(breaks = seq(0, .6, by=.3), labels = percent_format(accuracy = 2))  +
ylab("Accessible binding sites (IDR < 0.03)")

# Save PDF
ggsave("IMG/Fig4B_accessible_binding_sites.pdf", width = 3, height=1.25)
