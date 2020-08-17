#############################################################
#
# Cistrome datasets with highest overlap of the P53-C2 class
#
# Wout Megchelenbrink
# January 22, 2020
#
#############################################################

source("include/style.R")

DT
DT <- fread("DATA/Cistrome_H3K27ac_P53_C2.tsv.gz")[targets.pct >= 35]
setorder(DT, -JI)

DT[, treatment:="Unclassified"]
DT[GSMID %in%  c("GSM2438478", "GSM2438479", "GSM2197719", "GSM2197720", "GSM2197723"), treatment:="Nocodazole"]
DT[GSMID %in%  paste0("GSM16486", 52:55), treatment:="Reprogramming"]
DT[GSMID %in% paste0("GSM11044", 53:57), treatment:="Differentiation"]
DT[GSMID %in% c("GSM1355171","GSM1355172"), treatment:="Differentiation"] # Wysocka lab
DT[GSMID %in% paste0("GSM1410356"), treatment:="Differentiation"] # Krishnakumar CSC
DT[GSMID %in% c("GSM1618722", "GSM1618723"), treatment:="UV radiation"]

DT[GSMID %in% c("GSM1631257"), treatment:="Reprogramming"] # Jacob Hannah
DT[GSMID %in% paste0("GSM21984", 56:59), treatment:="Nutrient starvation"]
DT[GSMID %in% c("GSM1499128", "GSM1499129"), treatment:="Zic2 RNAi"]
DT[GSMID %in% c("GSM2417095"), treatment:="Reprogramming"] # iPSC Kathryn Plath
DT[GSMID %in% c("GSM2361227"), treatment:="Reprogramming"] # iPSC Thompson lab
DT[GSMID == "GSM2759401", treatment:="EOMES GFP reporter ESCs"]
DT[GSMID == "GSM1234533", treatment:="P300 inhibitor"]

DT <- DT[treatment != "Unclassified"]
DT[, GSMID:=factor(GSMID, levels=GSMID)]

DT[, cell.type:="Unclass"]
DT[Cell_type %in% c("Embryonic Stem Cell",  "Stem cell"),  cell.type:="ESC"]
DT[Cell_type == "Fibroblast", cell.type:="Fibroblast"]
DT[Cell_type == "Embryonic Fibroblast" | GSMID %in% c("GSM2198459", "GSM2198458"), cell.type:="MEF"]
DT[Cell_type %in% c("iPSC", "Intermediate") | GSMID %in% c("GSM1631257", "GSM2361227"), cell.type:="iPSC"]
DT[Cell_type == "Erythroid progenitor", cell.type:="Erythrocyte Prog."]

# Get the H3K27ac in 2iL/SL (baseline)
k27ac.gr <- makeGRangesFromDataFrame(fread("DATA/H3K27ac_2iL_SL_union_IDR_0.03.bed.gz"))
e1.gr <- makeGRangesFromDataFrame(fread("DATA/P53_ChIPseq_and_STARRseq.tsv.gz")[class_2iL == "C2" & class_SL == "C2", .(chr, start=motif_start-250, end=motif_end+250)])

ovl <- findOverlaps(e1.gr, k27ac.gr)

# Fraction of peaks with H3K27ac in 2iL/SL WT cells
baseline.x  <- length(unique(queryHits(ovl))) / length(e1.gr)
baseline.y  <- length(unique(subjectHits(ovl))) / length(k27ac.gr)

# Make the plot
DT[, GSMID:=factor(GSMID, levels=rev(GSMID))]
ggplot(DT, aes(x=hm.pct, y=targets.pct, color=treatment, shape=cell.type)) +
geom_point() +
coord_flip() +
scale_color_manual(values = c("#f442df", "red", colors.basic)) +
geom_hline(yintercept = baseline.x*100, cex = .2, linetype = "dashed", color = colors.basic[3]) +
geom_vline(xintercept = baseline.y*100, cex = .2, linetype = "dashed", color = colors.basic[3]) +
ylab("P53-C1 peaks overlapping a H3K27ac peak (%)") +
xlab("H3K27ac peaks overlapping a P53-C1 peak") +
theme_wout_minimal() +
scale_y_continuous(limits = c(20, 100), breaks = seq(20, 100, by=20)) +
ggtitle("H3K27ac peaks at P53-C1 binding sites")

# Save PDF
ggsave("IMG/Fig4F_P53_C2_overlapping_Cistrome_H3K27ac_peaks.svg", width = 5, height = 3)
