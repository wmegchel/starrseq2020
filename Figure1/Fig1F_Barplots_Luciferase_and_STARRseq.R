##################################################################
#
# Barplots of STARR-seq and luciferase for 4 selected loci
# Wout Megchelenbrink
#
# Januari, 20 2020
#
##################################################################
rm(list=ls())

library(ggthemes)
library(RColorBrewer)
show_col(colorblind_pal()(8))

source("include/style.R")
pc <- 0.5 # pseudo count

# Read luciferase data
DT <- fread("DATA/Luciferase_WT_2iL_SL.tsv.gz")

# Log2 scale STARR-seq
DT[, enrichment_2iL:=log2(enrichment_2iL_GM+pc)]
DT[, enrichment_SL:=log2(enrichment_SL_GM+pc)]

# Log2 scale luciferase
DT[, lucif_2iL:=log2((lucif_2i1*lucif_2i2*lucif_2i3)^(1/3)+pc)]
DT[, lucif_SL:=log2((lucif_SL1*lucif_SL2*lucif_SL3)^(1/3)+pc)]
DT[, batch:=factor(batch)]

# Linear regression after log2 transform, correcting for batch
fit.2iL <- lm(data = DT, formula = enrichment_2iL ~ lucif_2iL + batch)
DT[, `2iL1_BC`:=predict(fit.2iL, newdata = DT[, .(enrichment_2iL, lucif_2iL=log2(lucif_2i1+pc), batch)])]
DT[, `2iL2_BC`:=predict(fit.2iL, newdata = DT[, .(enrichment_2iL, lucif_2iL=log2(lucif_2i2+pc), batch)])]
DT[, `2iL3_BC`:=predict(fit.2iL, newdata = DT[, .(enrichment_2iL, lucif_2iL=log2(lucif_2i3+pc), batch)])]

fit.SL <- lm(data = DT, formula = enrichment_SL ~ lucif_SL + batch)
DT[, `SL1_BC`:=predict(fit.SL, newdata = DT[, .(enrichment_SL, lucif_SL=log2(lucif_SL1+pc), batch)])]
DT[, `SL2_BC`:=predict(fit.SL, newdata = DT[, .(enrichment_SL, lucif_SL=log2(lucif_SL2+pc), batch)])]
DT[, `SL3_BC`:=predict(fit.SL, newdata = DT[, .(enrichment_SL, lucif_SL=log2(lucif_SL3+pc), batch)])]

# Retransform to normal data
DT[, `2iL1_BC`:=2^`2iL1_BC`-pc]
DT[, `2iL2_BC`:=2^`2iL2_BC`-pc]
DT[, `2iL3_BC`:=2^`2iL3_BC`-pc]

DT[, `SL1_BC`:=2^`SL1_BC`-pc]
DT[, `SL2_BC`:=2^`SL2_BC`-pc]
DT[, `SL3_BC`:=2^`SL3_BC`-pc]

# Get mean and stdev for luciferase
DT[, luc.2i:=mean(c(`2iL1_BC`, `2iL2_BC`,`2iL3_BC`)), by=ID]
DT[, luc.2i.sd:=sd(c(`2iL1_BC`, `2iL2_BC`,`2iL3_BC`)), by=ID]
DT[, luc.ser.sd:=sd(c(`SL1_BC`, `SL2_BC`,`SL3_BC`)), by=ID]
DT[, luc.ser:=mean(c(`SL1_BC`, `SL2_BC`,`SL3_BC`)), by=ID]

# Get mean and stdev for STARR-seq
DT[, starr.2i:=mean(c(enrichment_2iL1, enrichment_2iL2)), by=ID]
DT[, starr.ser:=mean(c(enrichment_SL1, enrichment_SL2)), by=ID]
DT[, starr.2i.sd:=sd(c(enrichment_2iL1, enrichment_2iL2)), by=ID]
DT[, starr.ser.sd:=sd(c(enrichment_SL1, enrichment_SL2)), by=ID]

# Convert to long format
DTX <- melt.data.table(DT, id.vars = "ID", measure.vars = list(c("luc.2i", "luc.ser", "starr.2i", "starr.ser"),
                                                           c("luc.2i.sd", "luc.ser.sd", "luc.ser.sd", "starr.ser.sd")))

# Set variables
DTX[variable == 1, `:=` (condition="2iL", measure="Luciferase")]
DTX[variable == 2, `:=` (condition="SL", measure="Luciferase")]
DTX[variable == 3, `:=` (condition="2iL", measure="STARR")]
DTX[variable == 4, `:=` (condition="SL", measure="STARR")]
DTX[, variable:=NULL]
setnames(DTX, c("ID", "avg", "sd", "condition", "measure"))


DTX[, measure:=factor(measure, levels=c("STARR", "Luciferase"))]
DT[ID %in% c(5, 18, 42, 19) ]
# Create the plot
ggplot(DTX[ID %in% c(5, 18, 42, 19) ], aes(x=condition, y=avg, fill = measure)) +
geom_errorbar(aes(x=condition, ymin=avg, ymax=avg+sd), width = .2, position = position_dodge(.9)) +
geom_bar(stat = "identity", position = position_dodge(), color = "white") +
scale_fill_manual(values = colors.dark[c(5,4)]) +
xlab("Condition") +
ylab("Enrichment") +
scale_y_continuous(limits = c(0,20), breaks= seq(0,20, by=5)) +
theme_wout_minimal() +
theme(legend.title = element_blank(), legend.position = c(0.8, .8)) +
facet_wrap(~ID)

# Save PDF
ggsave("IMG/Fig1F_Barplots_Luciferase_and_STARRseq.pdf", width = 3, height = 4.5)
