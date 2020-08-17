####################################################################################
#
# Barplot of the P53-C1, P53-C2 classes and the percentage with STARR-seq
# enrichment (FC > 3) in 2iL and SL
#
# Wout Megchelenbrink
# January 22, 2020
#
####################################################################################

rm(list=ls())
source("include/style.R")

# Read the data
P53 <- fread("DATA/P53_ChIPseq_and_STARRseq.tsv.gz")

# Count guys with STARR-seq enrichment
DT <- rbind(P53[, .(.N, condition="2iL"), by=.(class=class_2iL, STARR=ifelse(padj_2iL1 < 0.05 & padj_2iL2 < 0.05 & enrichment_2iL >= 3, "Yes", "No"))],
            P53[, .(.N, condition="SL"), by=.(class=class_SL, STARR=ifelse(padj_SL1 < 0.05 & padj_SL2 < .05 & enrichment_SL >= 3, "Yes", "No"))])

DT[, class:=factor(class, levels = rev(c("C1", "C2")))]
DT[, STARR:=factor(STARR, levels = c("Yes", "No"))]
colors <- c(colors.light[5], colors.basic[5])

# Make the barplot
ggplot(DT, aes(x=class, y=N, fill = rev(STARR), label = N)) +
geom_bar(stat = "identity", color = "white", cex = .5) + 
coord_flip() + 
scale_fill_manual(values=colors, guide=F) +
scale_y_continuous(breaks = c(0, 1250, 2500)) +
geom_text() + 
facet_wrap(~condition, nrow = 2) +
theme_wout_minimal() 

# Save PDF
ggsave("IMG/Fig4B_Barplot_P53_C1_C2_active_inactive.pdf", width = 4, height = 1)

