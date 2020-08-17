#######################################################################################
#
#
# Lineplot of the GC-bias in STARR-seq input libraries compare to ChIP-seq libraries
#
# Wout Megchelenbrink
# January 23, 2020
#######################################################################################

source("include/style.R")
DT <- fread("DATA/GC_bias_STARRseq_and_ChIPseq.tsv.gz")

# Get total reads per sample
DT[, n.reads:=sum(READ_STARTS), by=.(sample, type)]

# Factorize
DT[, type:=factor(type, levels=c("STARR-seq input", "ChIP-seq input"))]
DT[, sample:=factor(sample, levels = c("2iL", "SL"))]

# Make the plot
ggplot(DT, aes(x=GC/100, y=READ_STARTS/n.reads, color = sample, linetype = type)) +
geom_line() +
xlab("GC% (per 100bp window)") +
ylab("Reads in bin (%)") +
scale_color_manual(values = colors.basic[c(1,3)]) +
scale_x_continuous(labels =  percent_format()) +
scale_y_continuous(labels =  percent_format()) +
theme_wout_minimal() +
theme(legend.position = c(0.8, 0.8))

# Save PDF
ggsave("IMG/FigS1D_STARRseq_GC_bias_vs_ChIPseq.pdf", width = 3.75, height = 3.5)
