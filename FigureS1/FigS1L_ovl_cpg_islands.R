#######################################################################################
#
#
# Compute overlap of STARR-seq C1, C2, C3 classes with annotated CpG islands
# Wout Megchelenbrink
#
# January 23, 2020
#######################################################################################
# Read data and convert to GRanges
starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")[, .(starr_peak_id, chr, start, end, class)]
starr[,.N, by=class]
cpg <- fread("DATA/CpG_islands_mm9.bed.gz")

starr.gr <- makeGRangesFromDataFrame(starr)
cpg.gr <- makeGRangesFromDataFrame(cpg)

# Overlap with CpG islands
ovl <- findOverlaps(starr.gr, cpg.gr)
DT <- cbind(starr[queryHits(ovl), ], cpg[subjectHits(ovl),])

# Compute fraction per class
DT <- merge(DT[, .N, by=class], starr[, .N, by=class], by="class", suffixes = c("_ovl", "_tot"))
DT[, frac:=N_ovl/N_tot]

# Make the plot
ggplot(DT, aes(x=class, y=frac, fill = class)) +
geom_bar(color = "white", stat="identity", show.legend = F) +
scale_fill_manual(values = colors.basic[c(2,1,3)]) + 
scale_y_continuous(labels = percent_format()) +
theme_wout_minimal() 

# Save PDF
ggsave("IMG/FigS1I_CpG_island_overlap_per_class.pdf", width = 2, height = 3)