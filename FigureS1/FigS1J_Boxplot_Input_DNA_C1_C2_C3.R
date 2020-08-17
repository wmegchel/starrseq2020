#######################################################################################
#
# Boxplots of the input DNA coverage at C1, C2 and C3 loci
#
# Wout Megchelenbrink
# Feb 14, 2020
#######################################################################################

# Read data
starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")
starr[, input_2iL:=log2(input_2iL1+input_2iL2+1)]
starr[, input_SL:=log2(input_SL1+input_SL2+1)]

# Conver to long
DT <- melt.data.table(starr, id.vars = c("starr_peak_id", "class"), measure.vars = c("input_2iL", "input_SL"), variable.name = "condition", value.name = "log2_reads")
DT[, condition:=str_replace(condition, "input_", "")]

# Make the plot
ggplot(DT, aes(x=class, y=log2_reads, color=class, fill=class)) +
geom_boxplot(position = position_dodge(), outlier.size = .25, notch=F) +
stat_summary(geom = "crossbar",position = position_dodge(), fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
scale_fill_manual(values = colors.basic[c(2,1,3)], guide=F) +
scale_color_manual(values = colors.basic[c(2,1,3)], guide=F) +
xlab("") +
ylab("Input reads per peak (log2)") +
facet_wrap(~condition) +
theme_wout_minimal() 

# Save PDF
ggsave("IMG/FigS1G_Boxplots_input_DNA_per_class.pdf", width = 3, height = 3)


# Stats
DT[, mean(log2_reads), by=.(condition,class)]  

classes <- c("C1", "C2", "C3")
for(cond in c("2iL", "SL"))
{
  for(i in 1:length(classes))
  {
    for(j in 1:length(classes))
    {
      if(j > i)
      {
        cat(sprintf("%s vs %s  p = %2.2E\n", classes[i], classes[j], 
              wilcox.test(DT[condition==cond & class==classes[i], log2_reads], 
                          DT[condition==cond & class==classes[j], log2_reads])$p.value))
        
      }
    }
  }
}  
  

