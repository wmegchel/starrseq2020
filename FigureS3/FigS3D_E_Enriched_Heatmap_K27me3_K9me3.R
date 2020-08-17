########################################################################
# 
# Enriched Heatmap of H3K27me3 and H3K9me3 in C1 and C2
# Wout Megchelenbrink
# 
# September 22, 2019
#
########################################################################
rm(list=ls())
gc()
source("include/style.R")

# Params
n.top <- Inf
frac.cutoff <- .9
mark <- "H3K9me3" # choose from {H3K9me3,  H3K27me3}

# Read peaks and heatmap data
peaks <- fread("DATA/STARRseq_OSN_RND_250bp.tsv.gz")
lst   <- readRDS("DATA/Heatmap_data_H3K27me3_H3K9me3.Rds")

# Sample subset (only for testing)
if(!is.infinite(n.top))
  peaks <- peaks[sample.int(nrow(peaks), size = n.top)]

marks <- c(sprintf("%s_2iL", mark),  sprintf("%s_SL", mark))
epi.marks <- marks

ehm.list <- NULL
classes <- peaks[, sort(unique(Class))]
partition <- peaks[, factor(Class, levels=classes)]
hm.clust <- Heatmap(partition, 
                    col = structure(colors.dark[c(2,1,3,6)], names = classes), 
                    name = "Cluster", show_row_names = F, width = unit(4, "mm"), show_heatmap_legend = F)

mats <- list()
list(marks)
i <- 1

for(epi.mark in marks)
{
  cat(sprintf("Preparing output for mark = %s\n", epi.mark))
  
  # Put matrix in same order as the annotation
  lst[[epi.mark]] <- merge(peaks[, .(peak_id=starr_peak_id)], lst[[epi.mark]], by="peak_id", sort=F)
  
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
    col.max <- max(round(quantile(rowMaxs(mats[[i]]), probs = frac.cutoff),1), round(quantile(rowMaxs(mats[[i-1]]), probs = frac.cutoff),1))
  else
    col.max <- max(round(quantile(rowMaxs(mats[[i]]), probs = frac.cutoff),1), round(quantile(rowMaxs(mats[[i+1]]), probs = frac.cutoff),1))

  ehm.list <- ehm.list + EnrichedHeatmap(mats[[i]], name = epi.mark,  use_raster = T, raster_quality = 3, raster_device = "CairoPNG", col = colorRamp2(c(0,  .5*col.max, col.max), c(hm.cols.dark[3], "white",  hm.cols.dark[5])), #use_raster = T , col = colorRamp2(c(0,  .75*col.max, col.max), c("white", colors[floor((i+1)/2)], "#012663")),
                                         column_title = epi.marks[i], na_col = "pink", row_order = order(enriched_score(mats[[i]])  , decreasing = TRUE), #- enriched_score(mats[[2]])
                                         column_title_gp = gpar(fontsize = 9, fontface = "bold"), width = unit(15, "mm"),  # - enriched_score(mats[[1]])
                                         axis_name = ifelse(i == 1, c("", "", ""), c("","","")),
                                         axis_name_gp = gpar(fontsize = 8), 
                                         top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = colors.basic[c(2, 1, 3,6)]), show_error=F, 
                                                                                                  yaxis = ifelse(i == length(marks), T, F), 
                                                                                                  axis_param = list(facing = "outside", side = "right"), 
                                                                                                  ylim = c(0, 1),
                                                                                                  height = unit(1.25, "cm"))),
                                         heatmap_legend_param = list(title = NULL, color_bar = "continuous", legend_direction = "horizontal",
                                                                     legend_width = unit(15, "mm"), grid_border = "#555555", at = c(0, col.max),
                                                                     labels = c("0", sprintf(">%2.1f", col.max))))
  
  i <- i+1
}

# width and height are in inches => 25.4 mm
pdf(file = sprintf("IMG/FigS3D_E_EnrichedHeatmap_%s.pdf", mark), width = 8, height = 7)
draw(ehm.list,  heatmap_legend_side = "bottom", padding = unit(c(0,0,0,0), "mm"),  split = partition) 
dev.off()