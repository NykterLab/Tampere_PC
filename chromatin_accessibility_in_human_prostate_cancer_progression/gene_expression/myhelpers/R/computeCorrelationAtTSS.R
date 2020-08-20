#' computeCorrelationAtTSS
#' @description  Compute Pearson and Spearman correlation coefficients for each gene given in atac_mat. Determines if gene is mirna or protein coding and compute correlation. If protein expression is available, compute correlation also between atac and protein levels.
#' @param atac_mat a matrix of atac signals at TSS. Gene id as rownames.
#' @param gene_mat gene expression matrix.
#' @param mirna_mat mirna expression matrix.
#' @param protein_mat protein expression matrix.
#' @param debug set debug mode. Compute correlation only for a small subset of genes.
#' @return matrix with gene id as rownames and 4 columns: pearson, spearman, protein_pearson, protein_spearman
#' @export
#' @examples
#'
#' @import parallel
computeCorrelationAtTSS <- function (atac_mat, 
                                     gene_mat, 
                                     mirna_mat, 
                                     protein_mat, 
                                     samples_blacklist = NULL, 
                                     debug = FALSE, 
                                     randomize = FALSE,
                                     # random_seed = 123,
                                     ncores = detectCores()) {


  worker <- function (i,
                      .atac_mat = atac_mat,
                      .gene_mat = gene_mat,
                      .mirna_mat = mirna_mat,
                      .proteins_mat = proteins_mat,
                      .valid_samples = valid_samples,
                      .randomize = randomize) {

    res <- setNames(rep(NA, 8), c("pearson", "pearson_pval",
                                  "spearman", "spearman_pval",
                                  "pearson_protein", "pearson_protein_pval",
                                  "spearman_protein", "spearman_protein_pval"))

    gid <- rownames(.atac_mat)[i]
    atac <- .atac_mat[i,]

    if (gid %in% rownames(.mirna_mat)) {
      gene <- .mirna_mat[gid,]
      samples_key <- "mirna"
    } else {
      gene <- .gene_mat[gid,]
      samples_key <- "rna"
    }

    atac_gene <- atac[match(.valid_samples[[samples_key]], names(atac), nomatch = 0)]
    gene <- gene[match(.valid_samples[[samples_key]], names(gene), nomatch = 0)]

    if (sd(gene) > 0 & sd(atac_gene) > 0){
      
      
      if (.randomize) {
        atac_gene <- sample(atac_gene)
        gene <- sample(gene)
      }
      
      pear <- cor.test(atac_gene, gene, method = "pearson", alternative = "two.sided")
      spear <- cor.test(atac_gene, gene, method = "spearman", alternative = "two.sided")
      res[c(1,2,3,4)] <- c(round(pear$estimate, 2), pear$p.value,
                           round(spear$estimate, 2), spear$p.value)
      
    }

    if (gid %in% rownames(.proteins_mat)) {
      protein <- .proteins_mat[gid,]

      atac_protein <- atac[match(.valid_samples[["proteins"]], names(atac), nomatch = 0)]
      protein <- protein[match(.valid_samples[["proteins"]], names(protein), nomatch = 0)]

      tryCatch({
        if (sd(atac_protein) > 0 & sd(protein, na.rm=TRUE) > 0) {
          
          if(.randomize) {
            atac_protein <- sample(atac_protein)
            protein <- sample(protein)
          }
          
          pear_prot <- cor.test(atac_protein, protein, method = "pearson", alternative = "two.sided")
          spear_prot <- cor.test(atac_protein, protein, method = "spearman", alternative = "two.sided")

          res[c(5,6,7,8)] <- c(round(pear_prot$estimate, 2), pear_prot$p.value,
                               round(spear_prot$estimate, 2), spear_prot$p.value)
        }

      }, error = function (err) {
        print(cbind(atac_protein, protein))
        stop(err)
      })
    }
    return(res)
  }

  ##########
  ## MAIN ##
  ##########
  
  # set.seed(random_seed)

  valid_samples <- list(rna = intersect(colnames(atac_mat), colnames(gene_mat)),
                        mirna = intersect(colnames(atac_mat), colnames(mirna_mat)),
                        proteins = intersect(colnames(atac_mat), colnames(proteins_mat)))

  if (!is.null(samples_blacklist)) {
    valid_samples <- lapply(valid_samples, function(each) each[!each %in% samples_blacklist])
  }

  axis <- setNames(seq(nrow(atac_mat)), rownames(atac_mat))

  if (debug) {
    axis <- axis[sample(length(axis), 100)]
  }

  corlist <- mclapply(axis, worker, mc.cores = ncores)

  cormat <- do.call(rbind, corlist)

  return(cormat)
}
