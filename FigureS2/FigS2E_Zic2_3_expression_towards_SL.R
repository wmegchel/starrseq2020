#######################################################################################
#
#
# Lineplot (MEAN + SEM) of the RNA expresison of ZIC2 and ZIC3 in 2iL, SL and timepoints
# in between
#
# Wout Megchelenbrink
# January 23, 2020
#
#######################################################################################
rm(list=ls())
source("include/style.R")

# Read the data
DT <- fread("DATA/RNAseq_timeseries_2iL_SL.tsv.gz")[gene_name %in% c("Zic2" ,"Zic3", "Tcf3")]

# Measure variablegs
mv <- list("BR1"=c("2i_1", "SLD1_1", "SLD3_1", "SLD7_1", "SL_1", "EpiLC_72h_1"), 
           "BR2"=c("2i_2", "SLD1_2", "SLD3_2", "SLD7_2", "SL_2", "EpiLC_72h_2"))

conditions <- data.table(TP=factor(1:6), condition=c("2i", "SLD1", "SLD3", "SLD7", "SL", "EpiLC_72h"))
DTX <- melt.data.table(DT, id.vars = c("ensembl_gene", "gene_name"), measure.vars = mv , variable.name = c("TP"))
DTX <- merge(DTX, conditions, by="TP")
DTX[, TP:=as.integer(TP)]

# Compute average and SEM
DTX[, avg:=(BR1+BR2)/2]
DTX[,stdev:=sqrt((BR1 - avg)^2 + (BR2-avg)^2)]
DTX[, sem:=stdev / sqrt(2)]
DTX[, condition:=factor(condition, levels=c("2i","SLD1", "SLD3", "SLD7", "SL", "EpiLC_72h"))]
DTX <- DTX[!is.na(condition)]

DTX[, gene_name:=factor(gene_name, levels=c("Zic2", "Zic3", "Tcf3"))]

# Make the plot
ggplot(DTX, aes(x=TP, y=avg, color = gene_name)) +
geom_line() +
geom_errorbar(aes(x=TP, ymin=avg-sem, ymax=avg+sem), width = .25) +
scale_color_manual(values = colors.dark[1:3]) +
scale_x_continuous(breaks=1:6) +
scale_y_continuous(labels=comma, limits =  c(0,4200), breaks = seq(1000,4000, by= 1000)) +
theme_wout_minimal() 

# Save PDF
ggsave("IMG/FigS3E_Zic2_Zic3_RNAseq_2iL_to_SL_timeseries.pdf", width = 5, height = 4)  
