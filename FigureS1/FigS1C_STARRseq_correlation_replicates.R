#######################################################################################
#
#
# Scatterplot showing hte correlation between STARR-seq replicates (log2 scale)
# Uses the great rasterized geom_point from Viktor Petukhov's ggrastr package
#
# Wout Megchelenbrink
# January 23, 2020
#######################################################################################

library(ggrastr)

source("include/style.R")

# Read the data
starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")[class %in% c("C1","C2")]

# Convert to long format
DT <- rbind(starr[, .(starr_peak_id, s1=log2(enrichment_2iL1 + 1), s2=log2(enrichment_2iL2+1), condition = "2iL")],
            starr[, .(starr_peak_id, s1=log2(enrichment_SL1 + 1), s2=log2(enrichment_SL2+1), condition = "SL")])

## Make the plot
ggplot(DT, aes(x=s1, y=s2)) + 
ggrastr::geom_point_rast(cex = .01, color = colors.basic[1]) + 
scale_x_continuous(limits = c(-1,6.5), breaks = seq(0,6, by=2)) +
scale_y_continuous(limits = c(-1,6.5), breaks = seq(0,6, by=2)) +
scale_color_identity(guide = "none") +
geom_abline(slope = 1, color=colors.basic[3], linetype = "dashed") +
stat_density2d(aes(fill=..level.., alpha=..level..), geom="polygon", show.legend = T) +
scale_fill_gradient(low=colors.basic[2], high=colors.basic[3]) +
facet_wrap(~condition)  + 
theme_wout_minimal() +
xlab("Enrichment rep. 1 (log2)") +
ylab("Enrichment rep. 2 (log2)")

# Save PDF
ggsave("IMG/FigS1C_STARR_correlation_replicates.pdf", width = 5, height = 2.5)  

# Get PCC
DT[, cor(s1, s2),  by=condition]
