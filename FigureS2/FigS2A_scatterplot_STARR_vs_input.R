#######################################################################################
#
#
# Scatterplot showing correlation between input and STARR-seq AT STARR-seq peaks.
#
# Wout Megchelenbrink
# January 23, 2020
#
#######################################################################################
rm(list=ls())

library(cowplot)
library(ggExtra)
library(ggpubr)
source("include/style.R")

starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")[class %in% c("C1", "C2")]
DT <- rbind(starr[, .(starr_peak_id, input_DNA=log2(input_2iL1 + input_2iL2 + 1), RNA=log2(rna_2iL1 + rna_2iL2 + 1), condition="2iL")],
            starr[, .(starr_peak_id, input_DNA=log2(input_SL1 + input_SL2 + 1), RNA=log2(rna_SL1 + rna_SL2 + 1), condition="SL")])
            
scatterplot.with.margins <- function(DT)
{
  sp <- ggplot(DT, aes(x=input_DNA, y=RNA)) +
        ggrastr::geom_point_rast(size = .1, alpha = 0.6, shape = 20) + 
        scale_x_continuous(breaks = c(0,3,6,9,12)) +
        scale_y_continuous(breaks = c(0,3,6,9,12)) +
        geom_smooth(method = "lm", color = colors.basic[1]) +
        border()   +
        theme_wout_minimal() 

  # Marginal density plot of x (top panel) and y (right panel)
  xplot <- ggdensity(DT, "input_DNA")
  yplot <- ggdensity(DT, "RNA") +
  rotate()
  
  # Cleaning the plots
  sp <- sp + rremove("legend")
  yplot <- yplot + clean_theme() + rremove("legend") 
  xplot <- xplot + clean_theme() + rremove("legend")
  
  # Arrange the plot using cowplot
  gg <- plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
            rel_widths = c(2, 1), rel_heights = c(1, 2))

  return(gg)  
}

# Save plot to PDF
ggsave(plot = scatterplot.with.margins(DT[condition == "2iL"]), file = "IMG/FigS2A_scatter_input_vs_RNA_2iL.pdf", width = 4, height = 3)
ggsave(plot = scatterplot.with.margins(DT[condition == "SL"]), file = "IMG/FigS2A_scatter_input_vs_RNA_SL.pdf",  width = 4, height = 3)

# Correlation stats
unique(DT[, cor.test(input_DNA, RNA),by=condition], by="condition")[, .(condition, PCC=estimate, p.value)]
