#######################################################################################
#
#
# Barplot (MEAN + STDEV) of luciferase in WT and ZIC3 KO cells at 3 loci
#
# Wout Megchelenbrink
# January 23, 2020
#
#######################################################################################
source("include/style.R")

DT <- fread("DATA/Luciferase_WT_ZIC3KO_2iL_SL.tsv.gz")
DT[, starr_peak_id:=factor(starr_peak_id)]

# Make the barplot
ggplot(DT, aes(x=starr_peak_id, y=AVG, fill = genotype)) +
geom_bar(stat ="identity", position = "dodge") +
geom_errorbar(aes(x=starr_peak_id, ymin=AVG-STDEV, ymax=AVG+STDEV), position = position_dodge(width =  1), cex = .25, width = .5) +
facet_wrap(~condition, nrow = 2) +
scale_fill_manual(values=c("#FF7F00","#BBBBBB")) +
theme_wout_minimal() +
xlab(NULL) +
ylab("Normalized luciferase intensity") +
coord_flip()

# Save PDF
ggsave("IMG/Fig2G_Barplot_Luciferase_WT_Zic3KO.pdf", width = 6, height = 4)
