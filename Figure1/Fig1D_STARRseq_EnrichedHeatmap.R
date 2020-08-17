########################################################################
# 
# Makes an enriched heatmap of STARR-seq and histone modifcations commonly
# associated with enhancers
#
# Wout Megchelenbrink
# 
# January 20, 2020
#
########################################################################

rm(list=ls())
source("include/style.R")

# Params
max.quant       <- 0.75 # cutoff the gradient bar at this quantile for better visualization
n.bins.smooth   <- 5    # number of bins to smooth using running median
n.top           <- Inf # Set to number to subsample (mainly for testing)
classes <- c("C1", "C2", "C3")

# Get STARR-seq and histone modification strength in 120 bins of 50bp surrounding peak summit
lst <- readRDS("DATA/STARRseq_and_histone_modifications_log2_RPKM_6kb_50bp_bins.Rds")

length(lst)
for(i in 1:length(lst))
{
  cat(sprintf("i=%d, name=%s, nrow=%d\n", i, names(lst)[i], nrow(lst[[i]])))
}

# Get clusters
peaks <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")
peaks[, .N, by=class]

# Sample subset (only for testing)
if(!is.infinite(n.top))
  peaks <- peaks[sample.int(nrow(peaks), size = n.top)]

marks <- c("STARRseq_2i", "STARRseq_ser",
           "ATAC_2i","ATAC_ser",
           "H3K4me1_2i", "H3K4me1_ser", 
           "P300_2i", "P300_ser",
           "H3K27ac_2i", "H3K27ac_ser")
        
epi.marks <- marks
ehm.list <- NULL

partition <- peaks[, factor(class, levels=classes)]
hm.clust <- Heatmap(partition, 
                    col = structure(colors.dark[1:3], names = classes), 
                    name = "Cluster", show_row_names = F, width = unit(4, "mm"), show_heatmap_legend = F)

mats <- list()

i <- 1
for(epi.mark in marks)
{
  cat(sprintf("Preparing output for mark = %s\n", epi.mark))
  
  # Put matrix in same order as the annotation
  lst[[epi.mark]] <- merge(peaks[, .(starr_peak_id)], lst[[epi.mark]], by.x="starr_peak_id", by.y="peak_id", sort=F)
  
  mats[[i]]  <- data.matrix(lst[[epi.mark]][,-1])
  if(!epi.mark %in% names(lst))
    stop(sprintf("Mark %s not found", epi.mark))
  
  mats[[i]] <- pmax(mats[[i]],0) 
  
  # Smoothing: running medians of 5 bins
  mats[[i]] <- t(apply(mats[[i]], 1, function(x) runmed(x, k=n.bins.smooth)))
  
  class(mats[[i]]) <- c("normalizedMatrix", "matrix")
  attr(mats[[i]], "upstream_index") <- 1:60
  attr(mats[[i]], "target_index") <- integer(0)
  attr(mats[[i]], "downstream_index") <- 61:120
  attr(mats[[i]], "extend") <- c(3000, 3000)
  attr(mats[[i]], "target_is_single_point") <- FALSE 
  
  i <- i+1
}


i <- 1
for(epi.mark in marks)
{
  # Keep the color scaling for 2iL and SL the same
  if(i %% 2 == 0)
    col.max <- max(round(quantile(rowMaxs(mats[[i]]), probs = max.quant),1), round(quantile(rowMaxs(mats[[i-1]]), probs = max.quant),1))
  else
    col.max <- max(round(quantile(rowMaxs(mats[[i]]), probs = max.quant),1), round(quantile(rowMaxs(mats[[i+1]]), probs = max.quant),1))
  
  # Make the enriched heatmap; save bin colors as PNG to drastically save storage and loading time
  ehm.list <- ehm.list + EnrichedHeatmap(mats[[i]], name = epi.mark,  use_raster = T, raster_quality = 3, raster_device = "CairoPNG", col = colorRamp2(c(0,  .5*col.max, col.max), c(hm.cols.dark[3], "white",  hm.cols.dark[5])), #use_raster = T , col = colorRamp2(c(0,  .75*col.max, col.max), c("white", colors[floor((i+1)/2)], "#012663")),
                                         column_title = epi.marks[i], na_col = "pink", row_order = order(enriched_score(mats[[1]]) , decreasing = TRUE),
                                         column_title_gp = gpar(fontsize = 9, fontface = "bold"), width = unit(15, "mm"),  
                                         axis_name = ifelse(i == 1, c("", "", ""), c("","","")),
                                         axis_name_gp = gpar(fontsize = 8), 
                                         top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = colors.dark[1:3]), show_error=F, 
                                                                                                  height = unit(1.25, "cm"),
                                                                                                  yaxis = ifelse(i == length(marks), T, F), 
                                                                                                  axis_param = list(facing = "outside", side="right"),
                                                                                                  ylim = c(0,4.5))),
                                         gap = unit(1, "mm"),
                                         heatmap_legend_param = list(title = NULL, color_bar = "continuous", legend_direction = "horizontal",
                                                                     legend_width = unit(15, "mm"), grid_border = "#555555", at = c(0, col.max),
                                                                     labels = c("0", sprintf(">%2.1f", col.max))))
  
  i <- i+1
}

# Attach cluster row annotation
ehm.list <- ehm.list + hm.clust

# Save the PDF
pdf(file = "IMG/Fig1D_STARRseq_EnrichedHeatmap.pdf", width = 14, height = 7)
draw(ehm.list,  heatmap_legend_side = "bottom", padding = unit(c(0,0,0,0), "mm"),  split = partition) 
dev.off()
