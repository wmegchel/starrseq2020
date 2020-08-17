####################################################################
#
# Scatterplot of NRF1 occupancy in WT and DNMT3 TKO cells
#
# Wout Megchelenbrink
# January 22, 2020
####################################################################

rm(list=ls())
source("include/style.R")

# Rename class to C2
DT <- fread("DATA/NRF1_ChIPseq_WT_and_Dnmt3KO.tsv.gz")
DT[, WT:=log2(NRF1_WT_1 + NRF1_WT_2 +1)]
DT[, TKO:=log2(NRF1_TKO_1 + NRF1_TKO_2 +1)]

DT[, isDE:=0]
DT[padj < 0.05, isDE:=1]

# annotate the 4 strongly differential peaks
DT[starr_peak_id %in% c("starr_00375", "starr_18696", "starr_24521", "starr_29782"), label:=starr_peak_id]
DT[starr_peak_id %in% c("starr_00375", "starr_18696", "starr_24521", "starr_29782"), isDE:=2]
DT[, isDE:=factor(isDE)]
DT[, .N, by=isDE]

# Make the scatterplot
ggplot(DT, aes(x=WT, y=TKO, color=isDE, shape = isDE, label=label)) +
geom_point(cex = .75) +
scale_color_manual(values = c("#DDDDDD",colors.basic[1], colors.basic[3])) +
geom_abline(slope = 1, color = "#BBBBBB", linetype = "dashed") +
theme_wout_minimal() +
geom_text_repel(min.segment.length = 0, nudge_x = 1) +
scale_x_continuous(limits = c(0,10), breaks=c(0,5,10)) +
scale_y_continuous(limits = c(0,10), breaks=c(0,5,10)) +
xlab("NRF1 WT (log2 reads)") +
ylab("NRF1 DNMT3 TKO (log2 reads)")

ggsave("IMG/Fig4B_bottom_scatterplot_NRF1_WT_TKO.pdf", width = 4.25, height = 3)


DTX <- rbind(DT[padj < 0.05, .(.N, isDE="Yes"), by=class],
              DT[padj >= 0.05 | is.na(padj), .(.N, isDE="No"), by=class])

# Make the barplot
ggplot(DTX, aes(x=class, y=N, fill =isDE)) +
geom_bar(stat = "identity") +
scale_fill_manual(values = c(colors.light[6], colors.basic[1])) +
theme_wout_minimal() +
coord_flip() +
xlab(NULL) +
scale_y_continuous(limits = c(0,600), breaks = c(0,200,400,600)) +
ylab("STARR-seq enhancers (#)")

# Save PDF
ggsave("IMG/Fig3B_top_NRF1_enhancers_DE_non_DE.pdf", width = 3, height = 1)