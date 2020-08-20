#' subsetDarsInTadCorrelationTable
#' @description Perform the subsetting of correlation tables by finding dars in tad.
#' @param tad_gr GRanges representing TADs
#' @param dars_gr GRanges represeting dars or peaks
#' @param genes_gr GRanges representing genes
#' @param pearson_cormat Matrix of Pearson correlation coefficients between genes and atac features
#' @param pearson_pvalue Matrix of p-values associated with Pearson correlation coefficients between genes and atac features
#' @param spearman_cormat Matrix of Spearman correlation coefficients between genes and atac features
#' @param spearman_pvalue Matrix of p-values associated with Spearman correlation coefficients between genes and atac features
#' @param pearson_protein Matrix of Pearson correlation coefficients between proteins and atac features
#' @param pearson_protein_pvalue Matrix of p-values associated with Pearson correlation coefficients between proteins and atac features
#' @param spearman_protein Matrix of Spearman correlation coefficients between proteins and atac features
#' @param spearman_protein_pvalue Matrix of p-values associated with Spearman correlation coefficients between proteins and atac features
#' @import parallel
#' @import GenomicRanges
#' @export
#' @return a data.frame
subsetDarsInTadCorrelationTable <- function (tad_gr, dars_gr, genes_gr,
                                             pearson_cormat, pearson_pvalue,
                                             spearman_cormat, spearman_pvalue,
                                             pearson_protein_cormat, pearson_protein_pvalue,
                                             spearman_protein_cormat, spearman_protein_pvalue,
                                             debug = FALSE, ncores = detectCores()){

  worker <- function(i,
                     .tad_gr = tad_gr,
                     .dars_gr = dars_gr,
                     .genes_gr = genes_gr,
                     .dars_in_tad = dars_in_tad,
                     .genes_in_tad = genes_in_tad,
                     .pearson_cormat = pearson_cormat,
                     .spearman_cormat = spearman_cormat,
                     .pearson_pvalue = pearson_pvalue,
                     .spearman_pvalue = spearman_pvalue,
                     .pearson_protein_cormat = pearson_protein_cormat,
                     .spearman_protein_cormat = spearman_protein_cormat,
                     .pearson_protein_pvalue = pearson_protein_pvalue,
                     .spearman_protein_pvalue = spearman_protein_pvalue) {
    tryCatch({
      tad_id <- as.character(.tad_gr[i])

      dar_idx_in_current_tad <- .dars_in_tad[.dars_in_tad[,1] == i, 2]
      gene_idx_in_current_tad <- .genes_in_tad[.genes_in_tad[,1] == i, 2]

      dar_id_in_current_tad <- names(.dars_gr)[dar_idx_in_current_tad]
      gene_id_in_current_tad <- names(.genes_gr)[gene_idx_in_current_tad]

      if (length(dar_id_in_current_tad) > 0 & length(gene_id_in_current_tad) > 0){

        cur <- expand.grid(dar_id=dar_id_in_current_tad,
                           gene_id=gene_id_in_current_tad,
                           stringsAsFactors = FALSE)
        cur[,c("pearson", "pearson_pvalue",
               "pearson_protein", "pearson_protein_pvalue",
               "spearman", "spearman_pvalue",
               "spearman_protein", "spearman_protein_pvalue")] <- NA

        cur$tad_id <- tad_id

        for (j in seq(nrow(cur))) {

          dar_id <- cur[j, "dar_id"]
          gene_id <- cur[j, "gene_id"]

          cur[j, "pearson"] <- .pearson_cormat[dar_id, gene_id]
          cur[j, "spearman"] <- .spearman_cormat[dar_id, gene_id]

          cur[j, "pearson_pvalue"] <- .pearson_pvalue[dar_id, gene_id]
          cur[j, "spearman_pvalue"] <- .spearman_pvalue[dar_id, gene_id]

          if (gene_id %in% colnames(.pearson_protein_cormat)){
            cur[j, "pearson_protein"] <- .pearson_protein_cormat[dar_id, gene_id]
            cur[j, "spearman_protein"] <- .spearman_protein_cormat[dar_id, gene_id]

            cur[j, "pearson_protein_pvalue"] <- .pearson_protein_pvalue[dar_id, gene_id]
            cur[j, "spearman_protein_pvalue"] <- .spearman_protein_pvalue[dar_id, gene_id]
          }
        }
      } else {
        cur <- data.frame()
      }

    return(cur)
  }, error = function (error) {
    # cat(tad_id, dar_id, gene_id ,"\n")
    print(i)

    # print("dars_in_tad")
    # print(head(.dars_in_tad))
    #
    # print("genes_in_tad")
    # print(head(.genes_in_tad))

    print(error)
  })
  }

  if(debug) {
    dars_gr <- subset(dars_gr, names(dars_gr) %in% rownames(pearson_cormat))
    tad_gr <- tad_gr[tad_gr %over% dars_gr]
  }

  dars_in_tad <- findOverlaps(tad_gr, dars_gr)
  dars_in_tad <- as.matrix(dars_in_tad)

  genes_in_tad <- findOverlaps(tad_gr, genes_gr)
  genes_in_tad <- as.matrix(genes_in_tad)

  res <- do.call(rbind, mclapply(seq_along(tad_gr), worker, mc.cores = ncores))

  return(res)
}
