#######################################################################################
#
#
# Boxplots of GC-distribution per class
# Wout Megchelenbrink
#
# January 23, 2020
#######################################################################################

# Read the data
DT <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")

# Make the plot
ggplot(DT, aes(x=class, y=GC_center_500bp, fill = class, color = class)) +
geom_boxplot(notch = T, outlier.size = .5, show.legend =F) +
scale_fill_manual(values = colors.basic[c(2,1,3)]) + 
scale_color_manual(values = colors.basic[c(2,1,3)]) + 
stat_summary(geom = "crossbar", width=0.75, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
scale_y_continuous(labels = percent_format()) +
theme_wout_minimal() 

# Save PDF
ggsave("IMG/FigS1H_GC_percentage_per_class.pdf", width = 2, height =3)
