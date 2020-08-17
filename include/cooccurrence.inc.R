##########################################################################################
#
# Computes motif co-occurrence and removes overlapping (user defined width) motifs
# @author:        Wout Megchelenbrink
# @last update:   February 21, 2019
#
#
# Fixed fractional overlap, fixed the computation of overlapping motifs
# plot co-occurrence heatmap
# plot overlap density for a list of motifs
##########################################################################################

# dt <- fread("../Rtools/cooccurrence/wout_test_set_1000kmer_summits.txt", col.names = c("m_chr", "m_start", "m_end", "kmer", "chr", "start", "end"), select = 1:7)
# dtx <- dt[, .(kmer, chr = m_chr, start=m_start - start, end=m_end - start)]

# binding.position.plot <- function(dtx, kmer.list)
# {
#   dtx <- dtx[kmer %in% kmer.list]
#   dtx.ir <- IRanges(start=dtx$start, end=dtx$end)
#   bins.ir <- IRanges(start=1:100, end=1:100)
#   
#   ovl <- findOverlaps(bins.ir, dtx.ir)
#   dtq <- data.table(bin=queryHits(ovl) , kmer=dtx[subjectHits(ovl), kmer])
#   
#   hist <- dtq[, .N, by=.(kmer, bin)]
#   
#   ggplot(hist, aes(x=bin, y=N, col=kmer)) +
#   geom_line() 
# }


get.cooccurrence <- function(dt, max.overlap.nt = 5, max.overlap.frac = .25)
{
  dt <- add.unique.peakids(dt)
  mat.occ <- get.occurrence.matrix(dt)
  mat.ovl <- get.overlap.matrix(dt, max.overlap.nt, max.overlap.frac)
  
  cooc <- compute.cooccurrence(mat.occ, mat.ovl)
  
  return(cooc)  
}


compute.cooccurrence <- function(mat.occ, mat.ovl)
{
  O.unc <- t(mat.occ) %*% mat.occ           # observed co-occurrences
  O.cor <- O.unc - mat.ovl                  # observed co-occurrences corrected for overlap   
  
  # Number of occurrence per motif
  occ <- Matrix::colSums(mat.occ)   
  
  # The expected number of occurrences between motif A and B discarding overlapping motifs is
  # occurences of A wihout overlapping B / (# peaks where A does not overlap B ) * # occurences of B without overlapping A / (# peaks where A does not overlap B) = 
  # (occ(A) - ovl_AB)  (#peaks - ovl_AB)  * (occ(B) - ovl_AB) / (#peaks - ovl_AB) = 
  # occ(A)*occ(B) - (occ(A)+occ(B))*ovl_AB + ovl_AB^2 / #peaks^2 - 2*#peaks*ovl_AB + ovl_AB^2
  #
  #
  # Let's define
  # A = occ(A)*occ(B) 
  # B =  (occ(A)+occ(B))*ovl_AB
  # D =  #peaks^2 - 2*#peaks*ovl_AB + ovl_AB^2
  #
  #### Example
  # A has 914 occurrences,  B has 1137 occurrences, they overlap 554 times and we have 7072 peaks in total
  # The expected overlap free co-occurence =   (914-554) / (7072-554) * (1137-554) / (7072-554) = 0.004940175 =
  # [914*1137 - (914+1137)*554 + 554^2] /  7072^2 - 2*7072*554 + 554^2
  #
  #
  # In matrix language, we have
  D <- diag(occ)
  rownames(D) <- names(occ)
  colnames(D) <- names(occ)
  
  n.peaks <- nrow(mat.occ)
  A <- outer(occ, occ)
  B <- D %*% mat.ovl + t(D %*% t(mat.ovl)) 
  C <- n.peaks^2 - 2*n.peaks * mat.ovl  + mat.ovl^2
  
  E.cor <- (A - B + mat.ovl^2) / C * (n.peaks - mat.ovl)
  
  # A["Nanog(Homeobox)","Pitx1(Homeobox)"]
  # B["Nanog(Homeobox)","Pitx1(Homeobox)"]
  # E.cor["Nanog(Homeobox)","Pitx1(Homeobox)"]
  # 
  if(any(E.cor < 0))
    warning("Negative expectation found!")
  
  # Some statistics
  chi_2 <- ((O.cor - E.cor)^2) / E.cor
  chisq <- data.table(melt(data.matrix(chi_2), varnames = c("m1", "m2"), value.name = "chi_2"))
  ovl   <- data.table(melt(data.matrix(mat.ovl), varnames = c("m1", "m2"), value.name = "ovl"))
  O.unc <- data.table(melt(data.matrix(O.unc), varnames = c("m1", "m2"), value.name = "n_total"))
  O.cor <- data.table(melt(data.matrix(O.cor), varnames = c("m1", "m2"), value.name = "observed"))
  E.cor <- data.table(melt(data.matrix(E.cor), varnames = c("m1", "m2"), value.name = "expected"))
  dt <- cbind(O.unc, ovl[, .(ovl)], O.cor[, .(observed)], E.cor[, .(expected)], chisq[, .(chi_2)])
  
  # Margin motif/kmer counts
  dt <- merge(dt,  data.table(m1=names(occ), m1_count=as.integer(occ)), by="m1")
  dt <- merge(dt,  data.table(m2=names(occ), m2_count=as.integer(occ)), by="m2")
  
  # Total nr of peaks
  dt[, n_peaks:=nrow(mat.occ)]
  
  # Enrichment, pval and padj
  dt[, log2_enrich:=log2(observed + 1) - log2(expected + 1)]
  dt[, pval:=pchisq(chi_2, df=1, lower.tail = F)]
  dt[, padj:=p.adjust(pval, "BH")]
  
  setcolorder(dt, c("m1", "m2", "n_total", "ovl", "observed", "expected", "chi_2", "m1_count", "m2_count", "n_peaks", "log2_enrich", "pval", "padj"))
  
  dt <- dt[m1 != m2 & expected > 0]
  setorder(dt, padj)
  
  return(dt)
}


add.unique.peakids <- function(dt)
{
  # assign unique peak id
  dt[, peak_id:=.GRP, by=.(chr,start,end)]
  dt <- unique(dt, by=c("peak_id", "motif"))
  
  return(dt)
}


# Count overlapping motifs
get.overlap.matrix <- function(dtx, max.overlap.nt = Inf, max.overlap.frac = 1)
{
  motifs <- dtx[, sort(unique(motif))] 
  
  if(is.infinite(max.overlap.nt) & max.overlap.frac >= 1)
  {
    ovl <- Matrix(0, nrow = length(motifs), ncol = length(motifs), sparse = T, dimnames = list(motifs, motifs))
  } else
  {
    gr <- makeGRangesFromDataFrame(dtx, seqnames.field = "chr", start.field = "m_start", end.field = "m_end", ignore.strand = T)
    
    # Get the overlapping motifs 
    ovl <- findOverlaps(gr, gr)
    overlaps <- pintersect(gr[queryHits(ovl)], gr[subjectHits(ovl)])     
    
    # Fraction of overlap relative to the shortest motif
    frac.ovl <- width(overlaps) / pmin(width(gr[queryHits(ovl)]), width(gr[subjectHits(ovl)]))
    
    res <- cbind(dtx[queryHits(ovl), .(pid1=peak_id, m1=motif, m1_start=m_start, m1_end=m_end)], 
                 dtx[subjectHits(ovl), .(pid2=peak_id, m2=motif, m2_start=m_start, m2_end=m_end)], 
                 length.A=width(gr[queryHits(ovl)]),
                 length.B=width(gr[subjectHits(ovl)]),
                 ovl.nt=width(overlaps),
                 ovl.frac = frac.ovl,
                 coef=1)
    
    res <- unique(res, by=c("pid1", "pid2","m1", "m2"))
    
    # Keep the ones that exceed the cutoffs
    res <- res[ovl.nt > max.overlap.nt  | ovl.frac > max.overlap.frac]
    
    # Cast to sparse matrix
    mat <- dcast.data.table(data=res, formula = m1 ~ m2, fun.aggregate = sum, fill = 0, value.var = "coef")
    ovl <- Matrix(data.matrix(mat[, 2:ncol(mat)]), sparse = T, dimnames = list(mat$m1, colnames(mat)[2:ncol(mat)]))
    
    ovl <- ovl[motifs, motifs]  
  }
  
  return(ovl)
}


get.occurrence.matrix <- function(dt)
{
  dt[, occ:=1]
  occ <-  Matrix(acast(dt, peak_id ~ motif, value.var = "occ", fill = 0), sparse = T)
  
  return(occ)
}



get.cooccurrence.heatmap <- function(mat, triag="lower", row.title = "M1", col.title = "M2") 
{
  
  color.range = colorRamp2(c(low, 0, 20, 40, 60), c(colors.dark[3], "white", colors.lighter[4], colors.dark[4], colors.dark[5]))
  
  hm <- Heatmap(mat, col = color.range, rect_gp = gpar(border=T, col="#EEEEEE", lwd=.3),
                cell_fun =   function (j, i, x, y, width, height, fill) 
                {
                  if(i >= j)
                  {
                    if(triag == "lower" | triag == "both")
                    {
                      grid.rect(x, y, width, height, gp = gpar(fill = fill, border=T,  col="#DDDDDD", lwd=.3))    
                    } else
                    {
                      grid.rect(x, y, width, height, gp = gpar(fill = fill, border=T,  col="white", lwd=.3))  
                    }
                  }else
                  {
                    if(triag == "lower")
                    {
                      grid.rect(x, y, width, height, gp = gpar(fill = fill, border=T,  col="white", lwd=.3))  
                      
                    } else
                    {
                      grid.rect(x, y, width, height, gp = gpar(fill = fill, border=T,  col="#DDDDDD", lwd=.3))
                    }
                  }
                }, 
                row_title = row.title, column_title = col.title, 
                column_title_side = ifelse(triag == "upper", "top" , "bottom"), 
                row_title_side = ifelse(triag == "upper", "right" , "left"),
                top_annotation_height = unit(30, "mm"), 
                row_names_side =  ifelse(triag == "upper", "right" , "left"),
                column_names_side =  ifelse(triag == "upper", "top" , "bottom"),
                cluster_rows = T, cluster_columns = T, name = "Enrichment\nscore", 
                column_title_gp = gpar(fontsize=13), 
                row_title_gp = gpar(fontsize=13),
                column_names_gp = gpar(fontsize=4), 
                row_names_gp = gpar(fontsize=4)
  )
  
  return(hm)
}




