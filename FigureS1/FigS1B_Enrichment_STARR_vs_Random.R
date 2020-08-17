#######################################################################################
#
#
# Shows the enrichment in the STARR-seq (p < 0.05 in both replicates) compared to 
# 100,000 randomly selected regions with equal GC content
# 
# Wout Megchelenbrink
# January 23, 2020
#######################################################################################

source("include/style.R")

# Read the data
DT <- fread("DATA/STARRseq_signif_and_RND_enrichment.tsv.gz")

DT[type == "RND GC-matched" & enrichment >= log2(3), .N/1e5*100, by=condition]


# Make the plot
ggplot(DT, aes(x = enrichment, color = type, fill = type)) +
geom_density() +
facet_wrap(~condition) +
scale_fill_manual(values = alpha(colors.light[1:2], .5)) +
scale_color_manual(values = colors.basic[1:2]) +
geom_vline(xintercept = log2(3), linetype = "dashed", color = "#555555") +
xlab("Enrichment (log2)") +
theme_wout_minimal() +
theme(legend.position = c(.2,.9), legend.title = element_blank())

# Save PDF
ggsave("IMG/FigS1B_STARRseq_significant_vs_RND_GC_matched.pdf", width = 6, height = 4)
