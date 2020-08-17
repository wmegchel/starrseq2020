####################################################################
#
# CpG methylation in C1 and C2 regions with and witout NRF1 motif
#
# Wout Megchelenbrink
# January 22, 2020
####################################################################
rm(list=ls())

source("include/style.R")

# Get STARRseq classes and NRF1 motif occurrences
starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")[, .(starr_peak_id, chr, start=summit-250L, end=summit+250L, class)][class %in% c("C1", "C2")]
NRF1  <- fread("DATA/NRF1_motif_occurrences.tsv.gz")[, .(chr,start,end, NRF1=1L)]

# Annotate the STARR-seq with NRF1 motif
starr <- merge(starr, NRF1, by=c("chr", "start", "end"), all.x=T)
starr[is.na(NRF1), NRF1:=0L]
starr[, .N, by=.(class, NRF1)]

# Make GRanges
starr.gr <- makeGRangesFromDataFrame(starr)

lst <- list()
for(cond in c("2iL", "SL"))
{
  # CpG methylation data in 2iL and SL ESCs
  CpG <- fread(sprintf("DATA/CpG_DNA_Methylation_q10_rmdup_counts_%s_merged.tsv.gz", cond))
  
  # Make Granges
  gr <- makeGRangesFromDataFrame(CpG, start.field = "position", end.field = "position")
  
  # Overlap C1 and C2 regions with CpG meth.
  ovl <- findOverlaps(starr.gr, gr)
  
  # Concat
  lst[[cond]] <- cbind(starr[queryHits(ovl), .(starr_peak_id, chr, start, end, class, NRF1)], CpG[subjectHits(ovl), .(nMeth, nTot, culture=cond) ])
}

DT <- rbind(lst[["2iL"]], lst[["SL"]])

# Compute the total number of methylated / unmethylated CpG's per peak
DT <- DT[ , lapply(.SD, sum), by=.(starr_peak_id, chr, start, end, class, culture, NRF1), .SDcols=c("nMeth", "nTot") ]
DT[, length(unique(starr_peak_id)), by=.(class, NRF1)]

# Pct meth
DT[, pct:=nMeth/nTot * 100]

DT[, hasNRF1:=NULL]
DT[, hasNRF1:=paste0(class,"_", ifelse(NRF1==1, "Yes", "No"))]
DT[, hasNRF1:=factor(hasNRF1, levels=c("C1_Yes", "C1_No", "C2_Yes", "C2_No"))]

wilcox.test(DT[hasNRF1 == "C2_No" & culture == "2iL", pct], DT[hasNRF1 == "C2_Yes" & culture == "2iL", pct])$p.value
wilcox.test(DT[hasNRF1 == "C2_No" & culture == "SL", pct], DT[hasNRF1 == "C2_Yes" & culture == "SL", pct])$p.value

# Make the plot
ggplot(DT, aes(x=class, y=pct, fill = hasNRF1, color = hasNRF1)) + 
geom_boxplot(notch = F, outlier.shape =  NA,  cex = .25) +
stat_summary(geom = "crossbar", position = position_dodge2(width=0.5), fatten=0, cex=.25, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
theme_wout_minimal() +
scale_fill_manual(values = c(colors.basic[2], colors.light[2], colors.basic[1], colors.light[1])) +
scale_color_manual(values = c(colors.basic[2], colors.light[2], colors.basic[1], colors.light[1])) +
xlab(NULL) +
ylab("CpG methylation (%)") +
facet_wrap(~culture)

# Save PDF
ggsave("IMG/Fig3A_right_DNA_meth_C1_C2_with_NRF1_motif.pdf", width = 5, height = 2)
