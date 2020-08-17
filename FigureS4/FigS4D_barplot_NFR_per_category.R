####################################################################
#
# Nucleosome free regions at the 500bp center of P53-C1, P53-C2, OSN and random regions
#
# Wout Megchelenbrink
# Jan 30, 2020
####################################################################
rm(list=ls())
source("include/style.R")

# Read peaks
peaks <- fread("DATA/P53_OSN_RND.saf.gz")[class != "SWITCH"]
peaks.gr <- makeGRangesFromDataFrame(peaks, seqnames.field = "chr", start.field =  "start_250", end.field = "end_250")

# Read nucleosome free regions
NFR <- fread("DATA/NFR.tsv.gz")
NFR.gr <- makeGRangesFromDataFrame(NFR)

# Overlap and concat
ovl <- findOverlaps(peaks.gr, NFR.gr, maxgap = 0)
DT <- cbind(peaks[queryHits(ovl), ], NFR[subjectHits(ovl)])
DT <- unique(DT, by=c("peak_id", "culture", "cutoff"))

# Compute fractions
DTX <- DT[, .N, by=.(culture, cutoff, class)]
DTX <- merge(DTX, peaks[, .N, by=class], by="class",suffixes=c(".nfr", ".total"))
DTX[, pct:=N.nfr/N.total]

# Make the barplot
ggplot(DTX, aes(x=class, y=pct, fill=class)) +
geom_bar(stat = "identity") +
facet_grid(culture ~ cutoff) +
scale_fill_manual(values = colors.basic[c(2,1,4,6)]) +
theme_wout_minimal()  +
scale_y_continuous(limits = c(0,1), breaks =c(0, .5, 1))

# Save PDF
ggsave("IMG/FigS4D_P53_OSN_RND_NFR.pdf", width = 4, height = 3)
