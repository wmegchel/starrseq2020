#######################################################################################
#
#
# ATAC vs STARR-seq and H3K27ac vs STARRseq at C1-loci
# Wout Megchelenbrink
#
# August 16, 2020
#######################################################################################

library(data.table)
library(ggplot2)
library(viridisLite)

rm(list=ls())
source("include/style.R")

examples <- c("starr_00268", "starr_45893",  "starr_18133", "starr_08288", "starr_17758",
              "starr_28828", "starr_33564", "starr_25427",  "starr_07417", "starr_09944")

# Read data
starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")[class == "C1"]
k27ac <- fread("STARR_and_RANDOM_2kb_histone.tsv.gz")[mark %in% c("H3K27ac_2i", "H3K27ac_ser")]
DT.k27ac <- dcast.data.table(k27ac, formula= GeneID + Chr + Start + End ~ mark, value.var = "log2_rpkm")
DT.k27ac <- merge(starr[, .(starr_peak_id, Chr=chr, Start=summit-2000L, End=summit+2000L, enrichment_2iL, enrichment_SL)], DT.k27ac, by=c("Chr", "Start", "End"))


DT.k27ac <- melt.data.table(DT.k27ac, id.vars = c("starr_peak_id", "Chr", "Start", "End"), 
                           measure.vars = list(c("enrichment_2iL", "enrichment_SL"), c("H3K27ac_2i", "H3K27ac_ser")),
                           variable.name = "condition", value.name = c("STARR", "K27ac"))

DT.k27ac[, condition:=ifelse(condition==1, "2iL", "SL")]
DT.k27ac[, cor(K27ac, log2(STARR+1)), by=condition]


# Plot
ggplot(DT.k27ac, aes(x=K27ac, y=log2(STARR+1))) +
ggrastr::geom_point_rast(cex = .01, color = colors.basic[6]) + 
scale_color_identity(guide = "none") +
stat_density2d(aes(fill=..level.., alpha=..level..), geom="polygon", show.legend = T) +
scale_fill_gradientn(colours = viridis(n=20)) +
geom_point(data = DT.k27ac[starr_peak_id %in% examples],aes(x=K27ac, y=log2(STARR+1)), color = colors.basic[3]) + 
facet_wrap(~condition) + 
theme_wout_minimal()

# Save
ggsave("FigS1H_C1_H3K27ac_vs_STARRseq.svg", width = 4.5, height = 3)
