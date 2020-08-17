#######################################################################################
#
#
# ATAC vs STARR-seq at C1-loci
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

# Read ATAC- and STARR-seq data
starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")[class == "C1"]
atac <- fread("STARR_and_RANDOM_250bp_TF_and_ATAC.tsv.gz")[mark %in% c("ATAC_2i", "ATAC_ser")]
DT.atac <- dcast.data.table(atac, formula= GeneID + Chr + Start + End ~ mark, value.var = "log2_rpkm")

# Merge
DT.atac <- merge(starr[, .(starr_peak_id, Chr=chr, Start=summit-250L, End=summit+250L, enrichment_2iL, enrichment_SL)], DT.atac, by=c("Chr", "Start", "End"))

# Convert to long format
DT.atac <- melt.data.table(DT.atac, id.vars = c("starr_peak_id", "Chr", "Start", "End"), 
                           measure.vars = list(c("enrichment_2iL", "enrichment_SL"), c("ATAC_2i", "ATAC_ser")),
                           variable.name = "condition", value.name = c("STARR", "ATAC"))

# Factorize
DT.atac[, condition:=ifelse(condition==1, "2iL", "SL")]

# Compute pearson R
DT.atac[, cor(ATAC, log2(STARR+1)), by=condition]


# Make scatterplot
ggplot(DT.atac, aes(x=ATAC, y=log2(STARR+1))) +
ggrastr::geom_point_rast(cex = .01, color = colors.basic[6]) + 
scale_color_identity(guide = "none") +
stat_density2d(aes(fill=..level.., alpha=..level..), geom="polygon", show.legend = T) +
scale_fill_gradientn(colours = viridis(n=20)) +
geom_point(data = DT.atac[starr_peak_id %in% examples],aes(x=ATAC, y=log2(STARR+1)), color = colors.basic[3]) + 
facet_wrap(~condition) + 
theme_wout_minimal()

ggsave("IMG/FigS1G_ATAC_vs_STARR.svg", width = 4.5, height = 3)
