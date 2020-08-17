##################################################################
#
# Homer2 motif enrichment in STARR-seq peaks
# Wout Megchelenbrink
#
# Januari, 20 2020
#
##################################################################
library(matrixStats)
library(stringr)
# Parameters
score.min       <- 50   # min log10 pval
score.max       <- 300  # max log10 pval
min.fc          <- 1.5
min.pct.targets <- 5
culture.condition <- "SL"



# Read Homer2 output
DT  <- fread("Fig1_Homer_motifs.tsv", col.names = c("motif", "consensus", "pval", "logp", "log2_enrichment",
                                                    "n_peaks_with_motif", "pct.targets",
                                                    "n_bg_with_motif", "pct_bg_with_motif",
                                                    "class", "condition"))

DT[, pct.targets:=as.numeric(str_sub(pct.targets, end=-2))]
DT[, log10_pval:=pmin(300,-log10(exp(logp)))]


#DTmotif in at least 10% of the peaks, enrichment 2-fold or more and p.adj < 1e-5
#motifs <- DT[pct.targets >= min.pct.targets & log2_enrichment >= log2(min.fc) & log10_pval >= score.min , motif]
#sort(motifs)

#DT <- DT[DT[, .I[log10_pval >= max(log10_pval)], by=.(motif, condition)]$V1,]
DT <- unique(DT, by=c("motif", "class", "condition"))

# Convert to wide format
DTX <- dcast.data.table(DT[condition == culture.condition], formula = motif ~ class, value.var = "log10_pval")


DTX[motif %in% c("Tcf4(HMG)/Hct116-Tcf4-ChIP-Seq(SRA012054)/Homer",
                 "Tcf3(HMG)/mES-Tcf3-ChIP-Seq(GSE11724)/Homer",
                 "LEF1(HMG)/H1-LEF1-ChIP-Seq(GSE64758)/Homer"), TF:="TCF/LEF"]
                 

# Cluster similar motifs
DTX[motif %in% c( "KLF14(Zf)/HEK293-KLF14.GFP-ChIP-Seq(GSE58341)/Homer",                    
                  "EKLF(Zf)/Erythrocyte-Klf1-ChIP-Seq(GSE20478)/Homer", 
                  "KLF5(Zf)/LoVo-KLF5-ChIP-Seq(GSE49402)/Homer",
                  "KLF3(Zf)/MEF-Klf3-ChIP-Seq(GSE44748)/Homer", 
                  "Klf4(Zf)/mES-Klf4-ChIP-Seq(GSE11431)/Homer" ,
                  "MafK(bZIP)/C2C12-MafK-ChIP-Seq(GSE36030)/Homer",
                  "KLF6(Zf)/PDAC-KLF6-ChIP-Seq(GSE64557)/Homer",                                                       
                  "Klf9(Zf)/GBM-Klf9-ChIP-Seq(GSE62211)/Homer",                                                      
                  "Maz(Zf)/HepG2-Maz-ChIP-Seq(GSE31477)/Homer" ,                             
                  "Sp1(Zf)/Promoter/Homer",
                  "Sp5(Zf)/mES-Sp5.Flag-ChIP-Seq(GSE72989)/Homer",  
                  "Sp2(Zf)/HEK293-Sp2.eGFP-ChIP-Seq(Encode)/Homer"), TF:="KLF/SP"]

DTX[motif %in% c("Bach1(bZIP)/K562-Bach1-ChIP-Seq(GSE31477)/Homer",
                 "Bach2(bZIP)/OCILy7-Bach2-ChIP-Seq(GSE44420)/Homer",
                 "Fosl2(bZIP)/3T3L1-Fosl2-ChIP-Seq(GSE56872)/Homer",
                 "Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer",
                 "Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer",
                 "MafB(bZIP)/BMM-Mafb-ChIP-Seq(GSE75722)/Homer",
                 "NF-E2(bZIP)/K562-NFE2-ChIP-Seq(GSE31477)/Homer",
                 "NFE2L2(bZIP)/HepG2-NFE2L2-ChIP-Seq(Encode)/Homer",
                 "Nrf2(bZIP)/Lymphoblast-Nrf2-ChIP-Seq(GSE37589)/Homer"), TF:="AP-1"]

DTX[motif %in% c("NRF1(NRF)/MCF7-NRF1-ChIP-Seq(Unpublished)/Homer", 
                 "NRF(NRF)/Promoter/Homer"), TF:="NRF1"]

DTX[motif %in% c("FOXP1(Forkhead)/H9-FOXP1-ChIP-Seq(GSE31006)/Homer", 
                 "Foxo3(Forkhead)/U2OS-Foxo3-ChIP-Seq(E-MTAB-2701)/Homer"), TF:="FOXP1/FOXO3"]

DTX[motif %in% c("Brn1(POU,Homeobox)/NPC-Brn1-ChIP-Seq(GSE35496)/Homer",
                 "Oct2(POU,Homeobox)/Bcell-Oct2-ChIP-Seq(GSE21512)/Homer",
                 "Oct4(POU,Homeobox)/mES-Oct4-ChIP-Seq(GSE11431)/Homer",
                 "Oct6(POU,Homeobox)/NPC-Pou3f1-ChIP-Seq(GSE35496)/Homer"), TF:="POU"]

DTX[motif == "Sox2(HMG)/mES-Sox2-ChIP-Seq(GSE11431)/Homer", TF:="SOX2"]

DTX[motif %in% c("TEAD(TEA)/Fibroblast-PU.1-ChIP-Seq(Unpublished)/Homer",
                 "TEAD1(TEAD)/HepG2-TEAD1-ChIP-Seq(Encode)/Homer",
                 "TEAD2(TEA)/Py2T-Tead2-ChIP-Seq(GSE55709)/Homer ",
                 "TEAD4(TEA)/Tropoblast-Tead4-ChIP-Seq(GSE37350)/Homer",
                 "TEAD2(TEA)/Py2T-Tead2-ChIP-Seq(GSE55709)/Homer"), TF:="TEAD"]

DTX[motif %in% c("RAR:RXR(NR),DR5/ES-RAR-ChIP-Seq(GSE56893)/Homer",
                 "RARg(NR)/ES-RARg-ChIP-Seq(GSE30538)/Homer",
                 "SF1(NR)/H295R-Nr5a1-ChIP-Seq(GSE44220)/Homer",
                 "Nr5a2(NR)/mES-Nr5a2-ChIP-Seq(GSE19019)/Homer",
                 "Esrrb(NR)/mES-Esrrb-ChIP-Seq(GSE11431)/Homer"), TF:="ESRRB/NR5A2"]

DTX[motif %in% c("E2F3(E2F)/MEF-E2F3-ChIP-Seq(GSE71376)/Homer"), TF:='E2F']
DTX[motif %in% c("ELF1(ETS)/Jurkat-ELF1-ChIP-Seq(SRA014231)/Homer",
                 "ETS(ETS)/Promoter/Homer",
                 "ETV4(ETS)/HepG2-ETV4-ChIP-Seq(ENCODE)/Homer",
                 "Elk1(ETS)/Hela-Elk1-ChIP-Seq(GSE31477)/Homer",
                 "Elk4(ETS)/Hela-Elk4-ChIP-Seq(GSE31477)/Homer"), TF:='ETS']

DTX[motif == "OCT4-SOX2-TCF-NANOG(POU,Homeobox,HMG)/mES-Oct4-ChIP-Seq(GSE11431)/Homer", TF:="OSN"]
DTX[motif == "Tbx20(T-box)/Heart-Tbx20-ChIP-Seq(GSE29636)/Homer", TF:="TBX"] 
DTX[motif == "Zfp281(Zf)/ES-Zfp281-ChIP-Seq(GSE81042)/Homer", TF:="ZFP281"]
DTX[motif == "Zic3(Zf)/mES-Zic3-ChIP-Seq(GSE37889)/Homer", TF:="ZIC2/3"] 
DTX[motif == "IRF3(IRF)/BMDM-Irf3-ChIP-Seq(GSE67343)/Homer", TF:="IRF"] 

DTX[motif %in% c("p53(p53)/mES-cMyc-ChIP-Seq(GSE11431)/Homer", 
                 "p53(p53)/Saos-p53-ChIP-Seq(GSE15780)/Homer", 
                 "p63(p53)/Keratinocyte-p63-ChIP-Seq(GSE17611)/Homer", 
                 "p73(p53)/Trachea-p73-ChIP-Seq(PRJNA310161)/Homer" ), TF:="P53/63/73"]
DTX[motif == "STAT1(Stat)/HelaS3-STAT1-ChIP-Seq(GSE12782)/Homer", TF:="STAT"]

DTX[motif == "YY1(Zf)/Promoter/Homer", TF:="YY1"]

DTX[motif == "Tcfcp2l1(CP2)/mES-Tcfcp2l1-ChIP-Seq(GSE11431)/Homer", TF:="TCFCP2L1"]


# For each group of TFs, take the motif that is most enriched over the 3 clusters on average
DTX <- DTX[!is.na(TF)]
DTX[, tot:=rowSums2(data.matrix(DTX[, 2:4]))]
DTX <- DTX[DTX[, .I[tot >= max(tot)], by=TF]$V1,]
DTX <- unique(DTX,by="TF")

mat <- data.matrix(DTX[, .(C1,C2,C3)])
rownames(mat) <- DTX$TF

# Limit the score to at 1e-300
mat[mat > score.max] <- score.max

# Set color scale
col <- colorRamp2(c(0, score.max/2, score.max), c("white", colors.basic[4], colors.basic[3]))


TFs <- c("OSN", "YY1", "POU", "AP-1", "ESRRB/NR5A2",  "TCF/LEF", "ETS", "KLF/SP", "ZFP281", "ZIC2/3",
            "P53/63/73", "TBX", "IRF", "FOXP1/FOXO3", "NRF1")

# Make the plot
pdf(file = sprintf("IMG/Fig1H_Homer_Motifs_C1_C2_C3_%s.pdf", culture.condition), width = 3, height = 3.75)                 
hm <- Heatmap(mat[TFs, ], name = "Enrichment: -log10(pval)", col = col, #rect_gp = chm.border.style.medium, 
              cluster_columns = F, row_title_side = "left",
              cluster_rows = F,
              show_row_dend = F, column_title_side = "top",
              heatmap_legend_param = list( at=c(0, score.max/2,score.max),  legend_width = unit(35, "mm"), 
                                           legend_direction = "horizontal", title_position = "topcenter"))
ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom")
dev.off()