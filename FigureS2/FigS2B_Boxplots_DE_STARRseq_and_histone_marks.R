#######################################################################################
#
#
# Boxplots of TF occupancy and histone modifications at STARR-seq peaks that 
# are differential in 2iL or SL
#
# Wout Megchelenbrink
# January 23, 2020
#
#######################################################################################
rm(list=ls())
source("include/style.R")

# Read the data and get DE peaks
starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")[class == "C1" & abs(log2fc_SL_vs_2iL) >= log2(2.5) & padj_SL_vs_2iL < 0.05]
starr[, dir:=ifelse(log2fc_SL_vs_2iL > 0, "SL", "2iL")]
starr[, .N, by=dir]
marks <- fread("DATA/STARRseq_peaks_narrow_and_broad_marks.tsv.gz")

DT  <- merge(marks, starr[, .(starr_peak_id, dir)], by="starr_peak_id")

marks <- c("ATAC", "p300", "H3K27ac", "Med1", "Pou5f1", "Sox2", "Nanog") #, "Nr5a2", "Esrrb", "Klf4")
DT <- DT[mark %in% marks]
DT[, mark:=factor(mark, levels=marks)]

# Make TF and histone marks boxplots
ggplot(DT, aes(x=mark, y=log2_rpkm, color = condition, fill = condition)) +
geom_boxplot(notch = T, outlier.size = .25) + 
stat_summary(geom = "crossbar", position = position_dodge(width = 1),  fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
scale_color_manual(values = colors.basic[1:2]) +  
scale_fill_manual(values = colors.basic[1:2]) +  
facet_wrap(~dir, nrow = 2, strip.position = "right") +
theme_wout_minimal() +
theme(axis.text.x = element_text(angle = 90)) 

for(d in c("2iL", "SL"))
{
  cat("\n")
  for(m in sort(unique(DT$mark)))
  {
    pval <- wilcox.test(DT[mark == m & condition == "2iL" & dir == d, log2_rpkm], 
                        DT[mark == m & condition == "SL" &  dir == d, log2_rpkm],
                        alternative = "two.sided")$p.value
    
    cat(sprintf("Mark=%s, direction=%s, pval=%2.2E\n", m, d, pval))
  }
}
 


# Save PDF
ggsave("IMG/FigS2B_right_boxplots_TF_and_histone_marks_DE_STARRseq_peaks.pdf", width = 8, height = 4.5)


# Make STARR-seq boxplots
DTX <- rbind(starr[, .(starr_peak_id, mark = "STARR-seq", enrichment=log2(enrichment_2iL+1), condition = "2iL", dir)],
             starr[, .(starr_peak_id, mark = "STARR-seq", enrichment=log2(enrichment_SL+1), condition = "SL", dir)])

ggplot(DTX, aes(x=mark, y=enrichment, color = condition, fill = condition)) +
geom_boxplot(notch = T, outlier.size = .25) + 
stat_summary(geom = "crossbar", position = position_dodge(width = 1),  fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
scale_color_manual(values = colors.basic[1:2]) +  
scale_fill_manual(values = colors.basic[1:2]) +  
facet_wrap(~dir, nrow = 2, strip.position = "right") +
theme_wout_minimal() +
theme(axis.text.x = element_text(angle = 90)) +
xlab(NULL) +
ylab("Enrichment (log2)")


# Save PDF
ggsave("IMG/S2B_left_boxplots_STARR_DE_STARR.pdf", width = 2.75, height = 4.5)


for(d in c("2iL", "SL"))
{
  cat("\n")
  for(m in sort(unique(DTX$mark)))
  {
    pval <- wilcox.test(DTX[mark == m & condition == "2iL" & dir == d,enrichment], 
                        DTX[mark == m & condition == "SL" &  dir == d, enrichment],
                        alternative = "two.sided")$p.value
    
    cat(sprintf("Mark=%s, direction=%s, pval=%2.2E\n", m, d, pval))
  }
}

DTX[, length(unique(starr_peak_id)), by=dir]

