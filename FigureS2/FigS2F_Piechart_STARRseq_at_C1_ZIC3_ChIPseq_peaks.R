###################################################################################
#
# Piechart of the Zic3 ChIP-seq peaks with positive STARR-seq signal at promoters
# and enhancers
#
# Wout Megchelenbrink
# January, 20 2020
#
###################################################################################

rm(list=ls())
source("include/style.R")

# Read data
zic3  <- fread("DATA/ZIC3_ChIP_and_STARRseq.tsv.gz")
starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")[class %in% c("C1", "C2")]

# Convert to Grangws
zic3.gr <- makeGRangesFromDataFrame(zic3)
starr.gr <- makeGRangesFromDataFrame(starr)

# overlap
ovl <- findOverlaps(zic3.gr, starr.gr)

# Annotate those that overlap a STARR peak
zic3[, hasSTARR:=0]
zic3[queryHits(ovl), hasSTARR:=1]

# Count numbers and fractions
DT <- zic3[, .N, by=.(hasSTARR, isPMT)]
DTX <- DT[, .(hasSTARR, N, pct=N/sum(N) * 100), by=.(isPMT)]

# Factorize
DTX[, hasSTARR:=factor(hasSTARR)]
DTX[, isPMT:=ifelse(isPMT==1, "Promoter", "Enhancers")]
DTX[, isPMT:=factor(isPMT,levels=c("Promoter", "Enhancers"))]

# Make the promoter and enhancer pie charts
ggplot(DTX, aes(x="", y=pct, fill=hasSTARR)) + 
geom_bar(width = 1, stat = "identity", color="white") +
coord_polar("y", start=0) +
scale_fill_manual(values = rev(c(colors.basic[5], colors.light[5]))) +
theme(axis.text.x=element_blank()) +
theme_void() +
facet_wrap(~isPMT)

# Save PDF
ggsave("IMG/FigS2F_Piechart_ZIC3_STARR_PMT_ENH.pdf", width = 5, height = 2)