#' computeCorrelationDarsInPromoters
#' @description Compute correlation between genes, mirna, proteins and dars overlapping promoter ranges.
#' @param proms_gid GRanges of promoters
#' @param dars_gr GRanges of ATAC-seq features
#' @param atac_mat data.table of atac levels. Need dar_id column
#' @param gene_mat matrix of gene expression
#' @param mirna_mat matrix of mirna expression
#' @param proteins_mat matrix of proteins expression
#' @param samples_blacklist character vector of blacklisted samples
#' @param debug if debug, work on a smaller subset of overlaps
#' @param randomize boolean, shuffle samples order before correlation?
#' @param id_column name of the id column (dar_id or peaks_id)
#' @param ncores number of cores to use
#' @return data.frame
#' @export
#' @examples
#' @import GenomicRanges
#' @import parallel
computeCorrelationDarsInPromoters <- function (proms_gid,
                                               dars_gr,
                                               atac_mat,
                                               gene_mat,
                                               mirna_mat,
                                               proteins_mat,
                                               samples_blacklist = NULL,
                                               debug = FALSE,
                                               randomize = FALSE,
                                               # random_seed = 123,
                                               id_column = "dar_id",
                                               ncores = detectCores()) {

  worker <- function (i,
                      .olaps = olaps,
                      .proms_gid = proms_gid,
                      .dars_gr = dars_gr,
                      .atac_mat = atac_mat,
                      .gene_mat = gene_mat,
                      .mirna_mat = mirna_mat,
                      .proteins_mat = proteins_mat,
                      .valid_samples = valid_samples,
                      .randomize = randomize) {

    didx <- .olaps[i,1]
    gidx <- .olaps[i,2]

    gid <- names(.proms_gid)[gidx]
    did <- names(.dars_gr)[didx]

    # atac <- setNames(as.numeric(.atac_mat[dar_id == did, -"dar_id"]), names(.atac_mat)[-which(names(.atac_mat) == "dar_id")])
    atac <- .atac_mat[rownames(atac_mat) == did, ]

    tryCatch({
      if (gid %in% rownames(.mirna_mat)) {
        gene <- .mirna_mat[gid, ]
        key <- "mirna"
      } else {
        gene <- .gene_mat[gid, ]
        key <- "rna"
      }
    }, error = function (err) {
      print(i)
      print(gidx)
      print(gid)
      stop(err)
    })

    atac_gene <- atac[match(.valid_samples[[key]], names(atac))]
    gene <- gene[match(.valid_samples[[key]], names(gene))]

    if (.randomize) {
      atac_gene <- sample(atac_gene)
      gene <- sample(gene)
    }

    row <- data.frame(id=did, gene_id=gid,
                      pearson=NA, pearson_pvalue=NA,
                      spearman=NA, spearman_pvalue=NA,
                      pearson_protein=NA, pearson_protein_pvalue=NA,
                      spearman_protein=NA, spearman_protein_pvalue=NA,
                      stringsAsFactors = FALSE)


    if (sd(gene, na.rm = TRUE) > 0 & sd(atac_gene, na.rm = TRUE) > 0){

      pear <- cor.test(atac_gene, gene, method = "pearson", alternative = "two.sided")
      spear <- cor.test(atac_gene, gene, method = "spearman", alternative = "two.sided")

      row[,c("pearson", "pearson_pvalue", "spearman", "spearman_pvalue")] <- c(round(pear$estimate,2), pear$p.value,
                                                                               round(spear$estimate, 2), spear$p.value)
    }

    if (gid %in% .proteins_mat) {

      protein <- .proteins_mat[gid,]

      atac_protein <- atac[match(.valid_samples[["proteins"]], names(atac))]
      protein <-protein[match(.valid_samples[["proteins"]], names(protein))]

      if(.randomize) {
        atac_protein <- sample(atac_protein)
        protein <- sample(protein)
      }

      if (sd(protein, na.rm = TRUE) > 0 & sd(atac_protein, na.rm = TRUE) > 0){
        pear_prot <- cor.test(atac_protein, protein, method = "pearson", alternative = "two.sided")
        spear_prot <- cor.test(atac_protein, protein, method = "spearman", alternative = "two.sided")

        row[,c("pearson_protein", "pearson_protein_pvalue",
               "spearman_protein", "spearman_protein_pvalue")] <- c(round(pear_prot$estimate, 2), pear_prot$p.value,
                                                                    round(spear_prot$estimat, 2), spear_prot$p.value)
      }
    }

    return(row)
  }

  ##########
  ## MAIN ##
  ##########

  # set.seed(random_seed)

  valid_samples <- list(rna=intersect(colnames(atac_mat), colnames(gene_mat)),
                        mirna=intersect(colnames(atac_mat), colnames(mirna_mat)),
                        proteins=intersect(colnames(atac_mat), colnames(proteins_mat)))

  if (!is.null(samples_blacklist)) {
    valid_samples <- lapply(valid_samples, function(each) each[!each %in% samples_blacklist])
  }


  olaps <- findOverlaps(dars_gr, proms_gid)
  # olaps <- cbind(queryHits(olaps), subjectHits(olaps))
  olaps <- as.matrix(olaps)

  if (debug) {
    olaps <- olaps[sample(nrow(olaps), 100), ]
  }

  cormat <- do.call(rbind, mclapply(seq(nrow(olaps)), worker, mc.cores = ncores))
  colnames(cormat)[1] <- id_column
  rownames(cormat) <- NULL

  return(cormat)
}
