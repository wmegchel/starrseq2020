#######################################################################################
#
#
# Barplot (MEAN + STDEV) of luciferase in WT and TCF3 KO cells 
#
# Wout Megchelenbrink
# January 23, 2020
#
#######################################################################################
source("include/style.R")

# Read the data
DT <- fread("DATA/Luciferase_WT_TCF3KO_2iL_SL.tsv.gz")
DT[, genotype:=factor(genotype, levels=c("WT", "Tcf3KO"))]
DT[, starr_peak_id:=factor(starr_peak_id)]

# Make the barplot
ggplot(DT, aes(x=starr_peak_id, y=AVG, fill = genotype)) +
geom_bar(stat ="identity", position = "dodge") +
geom_errorbar(aes(x=starr_peak_id, ymin=AVG-STDEV, ymax=AVG+STDEV), position = position_dodge(width =  1), cex = .25, width = .5) +
facet_wrap(~condition, scales = "free") +
scale_fill_manual(values=c("#E20E32","#BBBBBB")) +
theme_wout_minimal() +
xlab(NULL) +
ylab("Normalized luciferase intensity") +
coord_flip()

# Save PDF
ggsave("IMG/FigS2E_Barplot_Luciferase_Tcf3KO_2iL.pdf", width = 5.5, height = 2.5)
