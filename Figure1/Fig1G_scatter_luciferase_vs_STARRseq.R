##################################################################
#
# STARR-seq vs Luciferase scatterplot
# Wout Megchelenbrink
#
# Januari, 20 2020
#
##################################################################

rm(list=ls())
source("include/style.R")
pc <- 0.5 # Pseudo count

# Read luciferade data
DT <- fread("DATA/Luciferase_WT_2iL_SL.tsv.gz")

# Log2 scale STARR-seq
DT[, enrichment_2iL:=log2(enrichment_2iL_GM+pc)]
DT[, enrichment_SL:=log2(enrichment_SL_GM+pc)]

# Log2 scale luciferase
DT[, lucif_2iL:=log2((lucif_2i1*lucif_2i2*lucif_2i3)^(1/3)+pc)]
DT[, lucif_SL:=log2((lucif_SL1*lucif_SL2*lucif_SL3)^(1/3)+pc)]
DT[, batch:=factor(batch)]

# Linear regression, correcting for batch
fit.2i <- lm(data = DT, formula = enrichment_2iL ~ lucif_2iL + batch)
DT[, enrichment_2i_batch_corr:=predict(fit.2i)]

fit.ser <- lm(data = DT, formula = enrichment_SL ~ lucif_SL + batch)
DT[, enrichment_SL_batch_corr:=predict(fit.ser)]

# Convert to long format
DTX <- melt.data.table(DT, id.vars = c("peak_id", "chr", "start", "end", "class"), 
                          measure.vars = list(c("enrichment_2iL", "enrichment_SL"), 
                                              c("enrichment_2i_batch_corr", "enrichment_SL_batch_corr")), 
                       variable.name = "condition", value.name = c("STARR", "Luciferase"))
DTX[, condition:=ifelse(condition == 1, "2iL", "SL")]

# Get correlation stats
DTX[, cor.test(Luciferase, STARR, method = "pearson"), by=condition]

# Borderline case that fell out of C1 with the new stricter threshold. 
# We will assign it to C1 manually here
DTX[peak_id == "AK_11053", class:="C1"]

# Make the scatter plot
ggplot(DTX, aes(x=STARR, y=Luciferase, color = class, label=peak_id)) + 
geom_point() +
geom_smooth(method = "lm", color = "black", fill = "black", alpha = .1) +
scale_color_manual(values = c("black", colors.dark[1:3])) +
scale_x_continuous(limits = c(-1,6), breaks = seq(0,6, by=2)) +
scale_y_continuous(limits = c(-1,6), breaks = seq(0,6, by=2)) +
facet_wrap(~condition, nrow = 2) +
xlab("STARR enrichment (log2)") + 
ylab("Luciferase enrichment (log2)") +
theme_wout_minimal() +
theme(legend.title = element_blank(), legend.position = "bottom")

# Save PDF
ggsave("IMG/Fig1G_Scatterplot_Luciferase_vs_STARRseq.svg", width = 3, height = 4.5)
