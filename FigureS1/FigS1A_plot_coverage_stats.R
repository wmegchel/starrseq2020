#######################################################################################
#
#
# STARR-seq genome wide input DNA coverage per sample
# 
# Wout Megchelenbrink
# January 23, 2020
#######################################################################################
rm(list=ls())
source("include/style.R")

# Cutoff coverage at 50 reads to prevent long and gapped tail
cov.max <- 50

# Read the  coverage data computed with  "bedtools genome cov"
DT <- fread("DATA/STARRseq_genome_cov.tsv.gz")
DT[cov > cov.max, cov:=cov.max]


DT[type == "input_DNA", sum(cov*frac)/4]

cov.mean <- DT[, .(cov.mean=sum(cov*frac)), by=.(sample, type)]

# Get mean coverage per sample/type
DTX <- merge(DT, cov.mean, by=c("sample", "type"))
DTX[cov > 0, cov.mean:=NA]
DTX[!is.na(cov.mean), label:=sprintf("mean = %2.1f", cov.mean)]

# Make the plot (only for input)
ggplot(DTX[type == "input_DNA"], aes(x=cov, y=frac)) + 
geom_bar(stat="identity", fill = colors.basic[1], color = "white") +
geom_vline(aes(xintercept = cov.mean), linetype = "dashed", color = colors.basic[3]) +
geom_text(aes(x=cov.mean, y=0.2, label=label), nudge_x = 20 ) +
facet_wrap(~sample) + 
xlab("Coverage depth") + 
ylab("Genome fraction") +
theme_wout_minimal()

# Save PDF
ggsave("IMG/FigS1_STARRseq_Input_Coverage.pdf", width=5, height=4)