################################################################################
#
# Conpute and show motif co-occurrence in the
# C2 STARR-seq class
#
# Wout Megchelenbrink
# January 22, 2020
#
################################################################################
rm(list=ls())
library(Matrix)
source("include/style.R")
source("include/cooccurrence.inc.R")

# Params
score.max   <- 150
score.min   <- -25 # negative means statistcally depleted, i.e. observed co-occurence < expected co-occurence

# chr, start, end, motif, m_start,m_end
fg   <- fread("DATA/C1_C2_C3_zoops.tsv.gz")

fg[, motif:=str_sub(motif, end=str_locate(motif, "/")[,1]-1)]
setnames(fg, "motif_start", "m_start")
setnames(fg, "motif_end", "m_end")

fg <- fg[class == "C2"]
fg <- fg[!is.na(motif)]

# Get cooccurrence; allow no overlap between motifs
cooc <- get.cooccurrence(fg, max.overlap.nt = 0, max.overlap.frac = 0)

# Module C1: focus on the motifs enirched in this class
c1.enriched <- c("Klf9(Zf)",
                 "CRE(bZIP)",
                 "p53(p53)", 
                 "Tbx20(T-box)", 
                 "IRF2(IRF)",
                 "FOXP1(Forkhead)",
                 "NRF1(NRF)",
                 "HNF4a(NR),DR1",
                 "IRF3(IRF)",
                 "Fra1(bZIP)")

# Convert to matrix
C1 <- cooc[m1 %in% c1.enriched & m2 %in% c1.enriched, ]
C1[, score:=pmin(score.max,-log10(padj)) * sign(log2_enrich)]
mat <- acast(C1, formula = m1 ~ m2, value.var = "score", fill = 0)

# Make the heatmap
col <- colorRamp2(c(score.min, 0, score.max), c(colors.basic[1], "white", colors.basic[3]))

pdf("IMG/Fig3E_C2_cooccurence.pdf", width = 4, height = 4)
Heatmap(mat, name = "-log10(pval)", col = col, 
        heatmap_legend_param = list( at=c(score.min, 75, 0, score.max), legend_height = unit(50, "mm"),  legend_width = unit(100, "mm"), 
                                     legend_direction = "vertical", title_position = "topcenter"))
dev.off()