####################################################################################
#
# NucleoATAC occupancy profiles at P53-C1, P53-C2,
# OSN and random regions flanked by 500 bp
#
# Wout Megchelenbrink
# January 22, 2020
#
####################################################################################

rm(list=ls())
gc()
flank <- 500

source("include/style.R")

get.nuc.occ.profile <- function(condition)
{
  # Read NucleoATAC data
  nucatac <- fread(sprintf("DATA/P53_OSN_1kb_flanked_%s.occ.bedgraph.gz", condition), col.names = c("chr", "start", "end", "occ"))
  
  # Read SAF, P53 ChIP-seq data and merge
  saf <- fread("DATA/P53_OSN_RND.saf.gz", col.names = c("peak_id", "chr", "start", "end"), select = 1:4)
  P53 <- fread("DATA/P53_ChIPseq_and_STARRseq.tsv.gz")[, .(peak_id, class=get(sprintf("class_%s", condition)), center=round(motif_start+motif_end)/2)]
  P53 <- merge(saf, P53, by="peak_id", all.x=T)
  
  # Take center; if not P53 (which takes the motif center) then take the peak center
  P53[is.na(center), center:=round((start+end)/2)]
  P53[is.na(class), class:=str_sub(peak_id, end=str_locate(peak_id, "_")[,1]-1)]
  
  # Convert to GRanges
  nucactac.gr <- makeGRangesFromDataFrame(nucatac)
  P53.gr <- makeGRangesFromDataFrame(P53, start.field = "center", end.field = "center")
  
  # Overlap center with 500bp flank
  ovl <- findOverlaps(P53.gr, nucactac.gr, maxgap = flank)
  
  # Merge
  DT <- cbind(P53[queryHits(ovl), .(peak_id, class, center)], nucatac[subjectHits(ovl), ])[end - start == 1]
  DT[, dist:=start - center]
  DT <- DT[abs(dist) <= flank]
  
  DTX <- DT[, .(avg=mean(occ)), by=.(class, dist)]
  RND <- DT[class == "RND", .(avg.rnd=mean(occ)), by=.(dist)]
  
  DTY <- merge(DTX, RND, by="dist")
  DTY[, avg:=avg-avg.rnd]
  
  setorder(DTY, class, dist)
  
  # Take running median of 11 nt (that have data)
  DTY[, avg.smooth:=runmed(avg, k=11), by=class]
  
  DTY[, condition:=condition]
  return(DTY)
}

# Get average nucleosome occupancy track per condition
DTX <- rbindlist(lapply(c("2iL", "SL"), function(condition) get.nuc.occ.profile(condition)))

# Factorize 
DTX <- DTX[class != "RND"] # used to normalize
DTX[, class:=factor(class, levels=c("C1","C2", "OSN"))]

# Make the plot
ggplot(DTX, aes(x=dist, y=avg.smooth, color=class)) + 
geom_line() + 
geom_vline(xintercept = 0, color = colors.basic[3], linetype = "dashed", cex=.2) +
geom_vline(xintercept = -120, color = colors.basic[6], linetype = "dashed", cex=.2) +
geom_vline(xintercept = 120, color = colors.basic[6], linetype = "dashed", cex=.2) +
theme_wout_minimal() +
scale_x_continuous(breaks = c(-500,0,500)) +
scale_color_manual(values = colors.basic[c(2:1, 4)]) +
xlab("Distance from binding site (bp)") +
ylab("Average nucleosome occupancy") +
facet_wrap(~condition, nrow = 1)  +
theme(axis.text.x = element_text(angle = 0, hjust = .5), legend.title = element_blank()) 

# Save PDF
ggsave("IMG/Fig4D_NucleoATAC_500bp_window_scaled_to_rnd.pdf", width = 4, height = 3)