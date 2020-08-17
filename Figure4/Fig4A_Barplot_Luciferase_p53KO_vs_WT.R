#######################################################################################
#
#
# Barplot (MEAN + STDEV) of luciferase in WT and P53 KO cells 
#
# Wout Megchelenbrink
# January 23, 2020
#
#######################################################################################
source("include/style.R")

DT <- fread("DATA/Luciferase_WT_P53KO_2iL_SL.tsv.gz")
DT[, genotype:=factor(genotype, levels=c("WT", "Trp53KO"))]

# Make the plot
ggplot(DT, aes(x=starr_peak_id, y=avg,  fill = genotype)) +
geom_bar(stat = "identity",  position = position_dodge(width = 1),  color = "white") +
geom_errorbar(aes(x=starr_peak_id, ymin=avg-sem, ymax=avg+sem, group = genotype), position = position_dodge(width = 1), width = .5, cex = .25) +
scale_fill_manual(values = c(colors.basic[4], "#CCCCCC")) +
facet_wrap(~culture, nrow = 2) + 
coord_flip() + 
theme_wout_minimal() +
theme(legend.position = "top") +
ylab("Firefly vs Renilla") +
xlab(NULL)

# Save PDF
ggsave("IMG/Fig4A_Boxplot_Luciferase_Trp53_KO_vs_WT.pdf", width = 3, height = 3.5)

# Paired t-test
t.test(DT[genotype == "WT" & culture == "2iL", avg], DT[genotype != "WT" & culture == "2iL", avg], paired=T)
t.test(DT[genotype == "WT" & culture == "SL", avg], DT[genotype != "WT" & culture == "SL", avg], paired=T)


