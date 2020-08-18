#' mergeWithEnsembl
#' @description
#' @param  cormat target data.table
#' @param ann.gene Ensembl annotation data.frame
#' @param mirna_map mapping between mirBase id and Ensembl gene id
#' @return merged data.table
#' @export
mergeWithEnsembl <- function (cormat, ann.gene, mirna_map, target_column = NULL){
  if (!is.null(target_column)) {
    setnames(cormat, target_column, "gene_id", skip_absent = TRUE)
  }

  mir <- merge(cormat[grepl("hsa", gene_id)], mirna_map, by.x = "gene_id", by.y = "mirna_name", all.x = TRUE)
  mir <- merge(mir, ann.gene, by.x = "mirna_id", by.y = "gene_id", all.x = TRUE)
  mir$mirna_id <- NULL

  tmp <- merge(cormat[!grepl("hsa", gene_id)], ann.gene, by.x = "gene_id", by.y="gene_id", all.x = TRUE)
  res <- rbind(mir, tmp)
  
  if (!is.null(target_column)) {
    setnames(cormat, "gene_id", target_column, skip_absent = TRUE)
  }
  

  return(res)
}
