#' removeOverlappingDarsBySpearmanCorrelationValue
#' @description Remove overlapping DARs correlating with the same gene by absolute Spearman correlation value
#' @param df complete correlation matrix 
#' @param comparison string representing one of the three comparisons in the analysis (pc_bph, crpc_pc, crpc_bph)
#' @return fitered cormat
#' @export
#' @examples
#'
#' @import data.table
#' @import GenomicRanges
removeOverlappingDarsByCorrelationValue <- function (df, comparison) {
  
  if (!is.data.table(df)) {
    df <- as.data.table(df)  
  }
  
  df$keep <- FALSE
  genes <- unique(df[comparison_group == comparison, gene_id])
  
  for ( i in seq_along(genes)) {
    g <- genes[i]
    row_idx <- which(df$gene_id == g & df$comparison_group == comparison)
    
    if (length(row_idx) > 1){
      sbst <- df[row_idx]
      dars <- data.frame(chrom=gsub("(chr.+)_[0-9]+_[0-9]+", "\\1", sbst$dar_id), 
                         start=gsub("chr.+_([0-9]+)_[0-9]+", "\\1", sbst$dar_id), 
                         end=gsub("chr.+_[0-9]+_([0-9]+)", "\\1", sbst$dar_id),
                         pearson=sbst$pearson, 
                         spearman=sbst$spearman,
                         original_row_index=row_idx)
      dars <- makeGRangesFromDataFrame(dars, 
                                       starts.in.df.are.0based = FALSE, 
                                       keep.extra.columns = TRUE)
      olaps <- findOverlaps(dars)
      olaps <- as.matrix(olaps)
      olaps <- olaps[olaps[,1] != olaps[,2],]
      
      if(nrow(olaps) > 0) {
        idx <- olaps[,1]
        overlapping_dars <- dars[idx]
        
        keep.pearson <- which.max(abs(overlapping_dars$pearson))
        keep.spearman <- which.max(abs(overlapping_dars$spearman))

        maxs <- c(overlapping_dars[keep.pearson]$pearson, 
                  overlapping_dars[keep.spearman]$spearman)

        k <- c(overlapping_dars$original_row_index[keep.pearson], 
               overlapping_dars$original_row_index[keep.spearman])

        row_idx <- k[which.max(maxs)]
        
        # row_idx <- row_idx[keep.spearman]
        
      } 
    } 
    
    df[row_idx, keep := TRUE]
  }
  
  ret <- df[keep == TRUE]
  ret[, keep := NULL]
  return(ret)
  
}
