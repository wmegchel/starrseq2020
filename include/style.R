################ COLORS ################
# colors.dark   <- c("#222222",  "#0F8F03", "#0B3FBA", "#FF7000", "#FB0E0E")
# colors.light  <- c("#555555", "#07B513", "#2062D4", "#FF9000", "#FC3728") 
# colors.lighter  <- c("#999999", "#3fbf47", "#3777e5", "#FFB000", "#FC3728")
# colors.lightest <- c("", "#D6FFD9")
# 
# hm.cols <- c("#e2371b", "#1b6cd6", "#61a01c", "#e8a51b", "#7248db", "#2c3138", "#2a8e75", "#2296c1", "#dd5b2c")
# 
# bs.col.dark <- c("#E9573F", "#4A89DC",  "#8CC152",  "#F6BB42",   "#967ADC", "#434A54", "#37BC9B",  "#3BAFDA",  "#DA4443", "#2062D4", "#FF9000")
# bs.col.light <- c("#FC6E51", "#5D9CEC", "#A0D468", "#FFCE54", "#AC92EC", "#656D78","#48CFAD",  "#4FC1E9", "#ED5565")
# bs.col.lighter <- c("#F48371", "#87BCF4", "#C6EA98", "#F9D991", "#b29ee2", "#6e7277")

################ GGPLOT 2 ################

#font_import()
#extrafont::loadfonts()

hm.cols <- c("#e2371b", "#1b6cd6", "#61a01c", "#e8a51b", "#7248db", "#2c3138", "#2a8e75", "#2296c1", "#dd5b2c")
hm.cols.dark   <- c("#222222",  "#0F8F03", "#0B3FBA", "#FF7000", "#FB0E0E")

colors.basic <- c(RColorBrewer::brewer.pal(n=10, "Paired")[seq(2,10, by=2)], "#555555")
colors.light <- c(RColorBrewer::brewer.pal(n=10, "Paired")[seq(1,10, by=2)], "#CCCCCC")
colors.dark <- c("#009E73", "#0072B2", "#D55E00", "#000000", "#CC79A7", "#F0E442", "#E69F00")


theme_wout  <- function (base_size = 10, base_family = "Arial", ...)
               {
                  theme_classic (base_size = base_size, base_family = base_family)  %+replace%
                    theme (
                      axis.line = element_line(size = .5, color = "black"),
                      axis.line.x = element_line(size = .5, color = "black"),
                      axis.line.y = element_line(size = .5, color = "black"),
                      panel.grid.major.x = element_line(size = 0.2, color = "#AAAAAA", linetype = "solid"),
                      panel.grid.minor.y = element_line(size = 0, color = NA), 
                      panel.grid.major.y = element_line(size = 0.2, color = "#AAAAAA", linetype = "solid"),
                      panel.grid.major = element_line(size = 0, color = NA, linetype = "solid"),
                      panel.grid.minor.x = element_line(size = 0, color = NA),
                      panel.background = element_rect(fill = NA, color = NA),
                      legend.background = element_rect(fill=NA, color=NA),
                      strip.text = element_text(size = base_size-1, colour = "black", face = "bold"), # facet wrap strip
                      strip.background = element_blank(), # element_rect(fill=NA, color=NA, size = 5), # no background in facet
                      panel.spacing = unit(0.5, "lines"), 
                      plot.title = element_text(size = 9),
                      
                #      legend.text=element_text(size=base_size-2),
                #      legend.title=element_text(size=base_size-1),
                  ) # little bit more distance between facets
               }


theme_wout_minimal  <- function (base_size = 10, base_family = "Helvetica", ...)
{
  theme_classic (base_size = base_size, base_family = base_family)  %+replace%
    theme (
      axis.line = element_line(size = .5, color = "black"),
      axis.line.x = element_line(size = .5, color = "black"),
      axis.line.y = element_line(size = .5, color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(), 
      panel.grid.major.y = element_blank(),
      panel.grid.major =element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.background = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill=NA, color=NA),
      strip.text = element_text(size = base_size-1, colour = "black", face = "bold"), # facet wrap strip
      strip.background = element_blank(), # element_rect(fill=NA, color=NA, size = 5), # no background in facet
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(size = 8, hjust = 0),
      plot.tag = element_text(face = "bold", size = 12)
    )
}


scientific_10 <- function(x) 
{
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

#?gpar
######## COMPLEX HEATMAP ##########
# hm <- list()
# hm$basesize <- 10
# #hm$width <- 50
# ht_global_opt(heatmap_row_names_gp = gpar(fontsize = hm$basesize), 
#               heatmap_column_names_gp = gpar(fontsize = hm$basesize),
#               heatmap_row_title_gp = gpar(fontsize = hm$basesize + 2, fontface = "bold"),
#               heatmap_column_title_gp = gpar(fontsize = hm$basesize + 2, fontface = "bold"),
#               heatmap_legend_labels_gp = gpar(fontsize = hm$basesize - 1),
#               heatmap_legend_title_gp = gpar(fontsize = hm$basesize - 1, fontface = "bold"))
# 
# border.style <- gpar(border=T, col="#DDDDDD", lwd=.2)
# chm.border.style.light <- gpar(border=T, col="#EEEEEE", lwd=.2)
# chm.border.style.medium <- gpar(border=T, col="#777777", lwd=.2)
# chm.border.style.dark <- gpar(border=T, col="#444444", lwd=.2)
# chm.border.style.white <- gpar(border=T, col="#FFFFFF", lwd=.2)
# chmm.dark <- c("#C60517", "#114F0E", "#0E2EC9", "#720970")
# chmm.light <- c("#D7505C", "#527D50", "#1284C4", "#9C529B")
# 

