########################################################################
# 
# Enriched Heatmap
# Wout Megchelenbrink
# 
# January 30, 2020
#
########################################################################

rm(list=ls())
source("include/style.R")

n.top <- Inf
classes <- c("C1", "C2", "OSN", "RND")

# Histone marks
lst <- readRDS("DATA/P53_Heatmap_data.Rds")
for(i in 1:length(lst))
{
  cat(sprintf("i=%d, name=%s, nrow=%d\n", i, names(lst)[i], nrow(lst[[i]])))
}

peaks <- fread("DATA/P53_OSN_RND.saf.gz")[class != "SWITCH"] # enhancers that switch from C1 to C2 or vice versa in 2iL <-> SL are omitted

# Sample subset (only for testing)
if(!is.infinite(n.top))
  peaks <- peaks[sample.int(nrow(peaks), size = n.top)]

marks <- c("STARRseq_2i", "STARRseq_ser",
           "P53_2i","P53_ser",
           "ATAC_2i","ATAC_ser",
           "DNAse1_2i","DNAse1_ser",
           "H3K4me1_2i", "H3K4me1_ser", 
           "H3K27ac_2i", "H3K27ac_ser",
           "P300_2i", "P300_ser",
           "Pol2_2i", "Pol2_ser")
epi.marks <- marks

partition <- peaks[, factor(class, levels=classes)]
hm.clust <- Heatmap(partition, 
                    col = structure(colors.basic[c(2,1,4,6)], names = classes), 
                    name = "Cluster", show_row_names = F, width = unit(4, "mm"), show_heatmap_legend = F)

colors <- viridis_pal()(length(marks))[c(1,2)]
colors <- rep(colors.basic[2], 10)

mats <- list()
i <- 1
for(epi.mark in marks)
{
  cat(sprintf("Preparing output for mark = %s\n", epi.mark))
  
  # Put matrix in same order as the annotation
  lst[[epi.mark]] <- merge(peaks[, .(peak_id)], lst[[epi.mark]], by="peak_id", sort=F)
  
  mats[[i]]  <- data.matrix(lst[[epi.mark]][,-1])
  if(!epi.mark %in% names(lst))
    stop(sprintf("Mark %s not found", epi.mark))
  
  mats[[i]] <- pmax(mats[[i]],0) 
  
  # Smoothing: running medians of 5 bins
  mats[[i]] <- t(apply(mats[[i]], 1, function(x) runmed(x, k=9)))
  class(mats[[i]]) <- c("normalizedMatrix", "matrix")
  attr(mats[[i]], "upstream_index") <- 1:60
  attr(mats[[i]], "target_index") <- integer(0)
  attr(mats[[i]], "downstream_index") <- 61:120
  attr(mats[[i]], "extend") <- c(3000, 3000)
  attr(mats[[i]], "target_is_single_point") <- FALSE 
  
  i <- i+1
}

ehm.list <- NULL
i <- 1
for(epi.mark in marks)
{
  # Set max color to 95%, 99% ? 
  if(i %% 2 == 0)
    col.max <- max(round(quantile(rowMaxs(mats[[i]]), probs = 0.95),1), round(quantile(rowMaxs(mats[[i-1]]), probs = 0.95),1))
  else
    col.max <- max(round(quantile(rowMaxs(mats[[i]]), probs = 0.95),1), round(quantile(rowMaxs(mats[[i+1]]), probs = 0.95),1))
  
  ehm.list <- ehm.list + EnrichedHeatmap(mats[[i]], name = epi.mark,  use_raster = T, raster_quality = 3, raster_device = "CairoPNG", col = colorRamp2(c(0,  .5*col.max, col.max), c(hm.cols.dark[3], "white",  hm.cols.dark[5])), #use_raster = T , col = colorRamp2(c(0,  .75*col.max, col.max), c("white", colors[floor((i+1)/2)], "#012663")),
                                         column_title = epi.marks[i], na_col = "pink", row_order = order(enriched_score(mats[[5]]) , decreasing = TRUE), #- enriched_score(mats[[2]])
                                         column_title_gp = gpar(fontsize = 9, fontface = "bold"), width = unit(15, "mm"),  
                                         axis_name = ifelse(i == 1, c("", "", ""), c("","","")),
                                         axis_name_gp = gpar(fontsize = 8), 
                                         top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = colors.basic[c(2,1,4,6)]), show_error=F, 
                                                                                                  yaxis = ifelse(i == length(marks), T, F), 
                                                                                                  axis_param = list(facing = "outside", side = "right"), 
                                                                                                  height = unit(1.25, "cm"),
                                                                                                  gap = unit(1, "mm"))),
                                         heatmap_legend_param = list(title = NULL, color_bar = "continuous", legend_direction = "horizontal",
                                                                     legend_width = unit(15, "mm"), grid_border = "#555555", at = c(0, col.max),
                                                                     labels = c("0", sprintf(">%2.1f", col.max))))
  
  i <- i+1
}

ehm.list <- ehm.list + hm.clust 

# width and height are in inches => 25.4 mm
pdf(file = "IMG/FigS4A_EnrichedHeatmap_P53_OSN_RND.pdf", width = 19, height = 7)
draw(ehm.list,  heatmap_legend_side = "bottom", padding = unit(c(0,0,0,0), "mm"),  split = partition) 
dev.off()