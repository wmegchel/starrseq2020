####################################################################
#
# Barplots showing the percentage of C1, C2, C2-motif, C3 and random
# regions that overlap a specifc type of ERV/retrotransposons
#
# Wout Megchelenbrink
#
# January 22, 2020
####################################################################

rm(list=ls())
source("include/style.R")

# Read the C1, C2, C3 classes (and the ones with the triple motifs)
DT <- fread("DATA/C1_C2_C3_RND_classes_ERVs.tsv.gz")[class != "C3"]
DT[,.N, by=class]

DT.gr <- makeGRangesFromDataFrame(DT)
class.sizes <- DT[,.N, by=class]

# Read the ERV regions
ERVs <- fread("DATA/ERVs_mm9_classified.tsv.gz")
ERVs.gr <- makeGRangesFromDataFrame(ERVs)
ovl <- findOverlaps(DT.gr, ERVs.gr)

# Get fraction per class/type that overlaps an ERV
DTX <- unique(cbind(DT[queryHits(ovl), ], ERVs[subjectHits(ovl), .(type) ]))
DTX <- DTX[, .N, by=.(type, class)]
DTX <- merge(DTX, class.sizes, by="class", suffixes = c(".ovl", ".total"))
DTX[, pct:=N.ovl/N.total * 100]

DTX[class=="C2_motif", class:="C2-CO"]
DTX[, class:=factor(class, levels=c("C2-CO", "C2", "C1", "RND"))]

# Make the barplots
ggplot(DTX, aes(x=type, y=pct)) +
geom_bar(stat="identity", fill=colors.basic[6]) +
facet_wrap(~class, ncol = 1, strip.position = "right") + 
scale_y_continuous(breaks = c(0,30,60)) +
theme_wout_minimal() +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save PDF
ggsave("IMG/Fig3F_STARRseq_ERV_barplots.pdf", width = 3, height = 3)  

# Statistics 
a <- DTX[class == "C2" & type=="IAP-I", N.ovl]
b <- DTX[class == "C2" & type=="IAP-I", N.total-N.ovl]
c <- DTX[class == "C2-CO" & type=="IAP-I",N.ovl]
d <- DTX[class == "C2-CO" & type=="IAP-I",N.total-N.ovl]

M <- matrix(c(a,b,c,d), nrow=2, byrow = F)
fisher.test(M)$p.value