#######################################################################################
#
# Subsampling of the merged STARR-seq libraries to estimate whether the number of enhancers
# increases much when the libraries would have been deeper sequenced
#
# Wout Megchelenbrink
# January 23, 2020
#######################################################################################

source("include/style.R")

# Read the data
DT <- fread("DATA/STARRseq_subsampling_peakcounts.tsv.gz")

# Convert to long format and compute mean + stdev
DTX <- melt.data.table(DT, id.vars = c("run", "frac"), value.name = "n.peaks", variable.name = "library")
DTX[library=="peak_2i", library:="2iL"]
DTX[library=="peaks_ser", library:="SL"]
DTX[library=="peaks_union", library:="union"]
DTX <- DTX[library %in% c("2iL", "SL")]
DTX <- DTX[, .(avg=mean(n.peaks), sd=sd(n.peaks)), by=.(frac, library)]

# Make the plot
ggplot(DTX[library %in% c("2iL", "SL")], aes(x=frac, y=avg, color = library)) +
geom_errorbar(aes(x=frac, ymin=avg-sd, ymax=avg+sd), width=.015) +
geom_line() +
theme_wout_minimal() +
xlab("Subsampling percentage of merged libraries") +
ylab("Significant peaks (#)") +
scale_color_manual(values=colors.basic[1:3])  +
geom_vline(xintercept = .5, linetype = "dashed") +
scale_x_continuous(labels = percent_format(), breaks = seq(.1,.9, by=.2)) +  
scale_y_continuous(labels = comma, limits=c(10000, 30000)) +
theme(legend.position = c(0.15, 0.8)) 

# Save PDF
ggsave("IMG/FigS1E_STARRseq_subsampling_saturation.pdf", width = 4, height = 3.5)
