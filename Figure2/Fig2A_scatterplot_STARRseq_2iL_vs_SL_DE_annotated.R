###################################################################################
#
# Scatterplot of STARR-seq signal in 2iL and SL, with DE peaks annotated
#
# Wout Megchelenbrink
# January, 20 2020
#
###################################################################################
rm(list=ls())
library(ggrastr)
source("include/style.R")

# Parameters
FC.cutoff <- 2.5 # fold change cutoff to consider a peak differential

# Read STARR-seq data
starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")[class %in% c("C1", "C2")]
starr[, isDE:=0]

# Annotate and color DE peaks
starr[padj_SL_vs_2iL < 0.05 & abs(log2fc_SL_vs_2iL) >= log2(FC.cutoff), isDE:=1]
starr[isDE == 0, col:="#CCCCCC"]
starr[isDE == 1 & log2fc_SL_vs_2iL > 0, col:="#9A9E04"]
starr[isDE == 1 & log2fc_SL_vs_2iL < 0, col:="#E20E32"]

### Make the barplot (top)
DT <- starr[class == "C1" & abs(log2fc_SL_vs_2iL) >= log2(FC.cutoff) & padj_SL_vs_2iL < .05, .N, by=.(condition=ifelse(sign(log2fc_SL_vs_2iL)==1, "SL", "2iL"))]
DT[, condition:=factor(condition, levels=c("SL", "2iL"))]
ggplot(DT, aes(x=condition, y=N, fill=condition)) +
geom_bar(stat = "identity") +
coord_flip() +  
xlab("") +
scale_fill_manual(values = c("#9A9E04", "#E20E32")) +
ylab("Differential peaks (#)") +
theme_wout_minimal()

ggsave("IMG/Fig2A_Top_Barplot_STARRseq_2iL_vs_SL_C1.pdf", width = 4, height = 1.5)  

#### Make the scatterplot (bottom)
ggplot(starr[class == "C1"], aes(x=log2(enrichment_2iL), y=log2(enrichment_SL), color = col)) +
geom_point_rast(cex = .01) + 
scale_x_continuous(limits = c(-1,6.5), breaks = seq(0,6, by=2)) +
scale_y_continuous(limits = c(-1,6.5), breaks = seq(0,6, by=2)) +
scale_color_identity(guide = "none") +
geom_abline(slope = 1, color = colors.dark[4], linetype = "dashed") +
xlab("STARR-seq 2iL (log2 enrichment)") +
ylab("STARR-seq SL (log2 enrichment)") +
theme_wout_minimal() +
facet_wrap(~class)

ggsave("IMG/Fig2A_Bottom_Scatterplot_STARRseq_2iL_vs_SL_C1.pdf", width = 4, height = 4)  
