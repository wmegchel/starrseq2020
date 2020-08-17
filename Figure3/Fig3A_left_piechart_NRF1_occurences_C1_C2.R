####################################################################
#
# Piechart of NRF1 motif occurrences at C1 and C2 peaks
# Wout Megchelenbrink
#
# January 22, 2020
####################################################################

# Get NR of STARR-seq peaks per class
starr <- fread("DATA/STARRseq_C1_C2_C3.tsv.gz")
class.sizes <- starr[, .N, by=class]

NRF1.occurences <- fread("DATA/NRF1_motif_occurrences.tsv.gz")
NRF1.occurences <- NRF1.occurences[, .N, by=class]

DT <- merge(class.sizes, NRF1.occurences, by  ="class", suffixes = c("_total", "_motif"))
DT[, N_no_motif:=N_total - N_motif]

# Convert to long format
DTY <- melt.data.table(DT, id.vars = "class", measure.vars = c("N_motif", "N_no_motif"))

DTY[, pct:=value/sum(value)*100, by=class]

# Make the plot
ggplot(DTY, aes(x="", y=pct, fill = variable)) +
geom_bar(width = 1, stat = 'identity', color = "white", cex = .25) +
coord_polar("y", direction = -1) +
facet_wrap(~class) +
theme_void() + 
scale_fill_manual(values = c(colors.basic[1], colors.light[1]))  

# Save PDF
ggsave("IMG/Fig3B_Piecahrt_NRF1_motif_occurrences_C1_C2.pdf", width = 3, height = 1.5)
