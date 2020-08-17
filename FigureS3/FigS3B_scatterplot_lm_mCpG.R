####################################################################
#
# Wout Megchelenbrink
# CpG methylation in C1, C2, C3 and random regions
# Figure S3B
#
# Jan 30, 2020
####################################################################

rm(list=ls())
source("include/style.R")
library(stringr)

# Homer motif presence per peak 
zoops <- fread("DATA/C1_C2_C3_zoops.tsv.gz")
zoops[, present:=1L]
zoops[, motif:=str_sub(motif, end=str_locate(motif, "/")[,1]-1)]
zoops <- unique(zoops, by=c("chr", "start", "end", "motif"))
zoops <- dcast.data.table(zoops, formula = chr + start + end + class ~ motif, value.var = "present", fill=0)

starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")[class %in% c("C1", "C2"), .(chr, start=summit-250L, end=summit+250L, class)]
starr.gr <- makeGRangesFromDataFrame(starr)

gc <- fread("DATA/CpG_DNA_Methylation_q10_rmdup_counts_SL_merged.tsv.gz")
gc[, start:=position]
gc[, end:=position+1]
gc.gr <- makeGRangesFromDataFrame(gc)

ovl <- findOverlaps(starr.gr, gc.gr)

DT <- cbind(starr[queryHits(ovl), ], gc[subjectHits(ovl), .(nMeth, nTot) ])
DT <- DT[, lapply(.SD, sum), by=.(chr, start, end), .SDcols=c("nMeth", "nTot")]
DT[, mCpG:=nMeth/nTot]

DT <- merge(DT[, .(chr, start, end, mCpG)], zoops, by=c("chr","start","end"))

# LM cannot handle names with brackets ..str_replace them
colnames(DT) <- str_replace_all(colnames(DT), "[(]", "_")
colnames(DT) <- str_replace_all(colnames(DT), "[)]", "_")
colnames(DT) <- str_replace_all(colnames(DT), "[,]", "_")
colnames(DT) <- str_replace_all(colnames(DT), "[.]", "_")

# Remove obsolete columns
DT[, chr:=NULL]
DT[, start:=NULL]
DT[, end:=NULL]
DT[, `NA`:=NULL]

mat <- matrix(nrow = ncol(DT), ncol = 3)
rownames(mat) <- colnames(DT)
colnames(mat) <- c("coef","pval", "pct_meth")

# Single var linear regression
for(i in 3:nrow(mat))
{
  mat[i ,1] <- summary(lm(formula = mCpG ~ get(colnames(DT)[i]), data = DT))$coef[2,1]
  mat[i ,2] <- summary(lm(formula = mCpG ~ get(colnames(DT)[i]), data = DT))$coef[2,4]
  mat[i, 3] <- DT[get(colnames(DT)[i]) == 1, mean(mCpG)]

  if(i %% 10 == 0)
    cat(sprintf("Done %d of %d\n", i, nrow(mat)))
}
DT
mat["Zic3_Zf_", ]
summary(lm(formula = mCpG ~ p53_p53_, data = DT))

# Prepare data and annotate some point
DTX <- as.data.table(mat[3:nrow(mat),], keep.rownames = "motif")
DTX[, pval:=pmax(pval, 1e-300)]
DTX[grep("p53", motif, ignore.case = T), TF:=motif]
DTX[grep("nrf1", motif, ignore.case = T), TF:=motif]
DTX[grep("klf", motif, ignore.case = T), TF:=motif]
DTX[grep("zic", motif, ignore.case = T), TF:=motif]

setorder(DTX, pval)
DTX[is.na(TF), col:="#CCCCCC"]
DTX[!is.na(TF), col:=colors.basic[1]]

# Make the plot
ggplot(DTX, aes(x=coef*100, y=-log10(pval), label=TF, color=col)) +
geom_point(cex = .2) +
geom_text_repel() +
#geom_vline(xintercept = DT[, mean(mCpG)], type = "dashed") +
scale_x_continuous(limits = c(-25, 25)) +
scale_y_continuous(limits = c(0, 310)) +
scale_color_identity() +
theme_wout_minimal()

# Save PDF 
ggsave("IMG/FigS3B_Scatter_motif_mCpG_regression.pdf", width = 3, height = 2.5)
