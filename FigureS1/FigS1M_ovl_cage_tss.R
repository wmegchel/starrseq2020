#######################################################################################
#
#
# Compute overlap of STARR-seq C1, C2, C3 classes with Phantom5 CAGE peaks that
# were classified as "TSS"
#
# Wout Megchelenbrink
# January 23, 2020
#
#######################################################################################
rm(list=ls())
source("include/style.R")

# Read STARR-seq
starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")[, .(starr_peak_id, chr, start=summit, end=summit, class)]
starr.gr <- makeGRangesFromDataFrame(starr)

# Read Phantom CAGE peaks and define promoter region
tss   <- fread("DATA/Phantom5_TSS_classified_CAGE_peaks_mm9.bed.gz")
tss[strand == "+", start:=start-2500L]
tss[strand == "+", end:=end+500L]
tss[strand == "-", start:=start-500L]
tss[strand == "-", end:=end+2500L]
tss.gr <- makeGRangesFromDataFrame(tss)

# Overlap
ovl.tss <- findOverlaps(starr.gr, tss.gr)

# Count overlap fraction
DT <- unique(starr[queryHits(ovl.tss), ], by="starr_peak_id")[, .N, by=class]
DT <- merge(DT, starr[,.N,by=class],by="class", suffixes = c("_tss", "_tot"))
DT[, frac:=N_tss/N_tot]

# Make the plot
ggplot(DT, aes(x=class, y=frac, fill = class)) +
geom_bar(stat = "identity") + 
scale_fill_manual(values = colors.basic[c(2,1,3)]) + 
#scale_y_continuous(labels = percent_format()) +
theme_wout_minimal() 

# Save PDF
ggsave("IMG/FigS1J_C1_C2_C3_ovl_Phantom_TSS_classified_CAGE_peak.pdf", width = 4, height = 3)  
