#' computeCorrelationDarInTad
#' @description  Compute Pearson and Spearman correlation coefficients for each gene-dar couple fallingi the same topologically associated domain (TAD). First perform overlap between dar GRanges object and genes GRanges object against TAD GRanges object. For each TAD start a thread, get dars and genes overlapping that TAD, perform cartesian product and for each generated couple compute Pearson and Spearman correlation coefficients. If gene id is detected in proteins data matrix, compute correlation coefficients with that as well.
#' @param tad_gr GRanges object representing TADs
#' @param dars_gr GRanges object representing DARs
#' @param genes_gr GRanges object representing genes
#' @param atac_mat Data.table object representing accessibility values (ATAC-seq)
#' @param gene_mat matrix representing gene expression values for protein coding genes
#' @param mirna_mat matrix representing mirna expression values
#' @param proteins_mat matrix representin protein expression values
#' @param debug Boolean. If TRUE perform operations on small subset of TAD (3)
#' @param ncores Number of available cores for multithread operations
#' @return matrix with gene id as rownames and 4 columns: pearson, spearman, protein_pearson, protein_spearman
#' @export
#' @examples
#'
#' @import parallel
#' @import GenomicRanges
computeCorrelationDarInTad <- function (tad_gr,
                                        dars_gr,
                                        genes_gr,
                                        atac_mat,
                                        gene_mat,
                                        mirna_mat,
                                        proteins_mat,
                                        debug = FALSE,
                                        randomize = FALSE,
                                        random.seed = 123,
                                        ncores = detectCores()) {

  worker <- function (i,
                      .dars_in_tad = dars_in_tad,
                      .genes_in_tad = genes_in_tad,
                      .tad_gr = tad_gr,
                      .gene_mat = gene_mat,
                      .mirna_mat = mirna_mat,
                      .proteins_mat = proteins_mat,
                      .atac_mat = atac_mat,
                      .valid_samples = valid_samples,
                      .randomize = randomize) {

    # find dar_id in current tad
    didx <- .dars_in_tad[.dars_in_tad[,2] == i, 1]
    dars_in_current_tad <- dars_gr[didx]
    dar_id_in_current_tad <- names(dars_in_current_tad)

    # find gene_id in current tad
    gidx <- .genes_in_tad[.genes_in_tad[,2] == i, 1]
    genes_in_current_tad <- genes_gr[gidx]
    gene_id_in_current_tad <- names(genes_in_current_tad)

    # do cartesian product between gene_id and dar_id
    current_cormat <- expand.grid(dar_id = dar_id_in_current_tad,
                                  gene_id = gene_id_in_current_tad,
                                  stringsAsFactors = FALSE)
    current_cormat[,"tad_id"] <- as.character(.tad_gr[i])
    current_cormat[,c("pearson", "spearman", "pearson_protein", "spearman_protein", "comparison")] <- NA

    # for each dar_id-gene_id couple compute correlation
    for (ii in seq(nrow(current_cormat))) {

      did <- current_cormat[ii,"dar_id"]
      gid <- current_cormat[ii, "gene_id"]

      tryCatch({
        if (gid %in% rownames(.mirna_mat)) {
          gene <- .mirna_mat[gid,]
          key <- "mirna"
        } else {
          gene <- .gene_mat[gid,]
          key <- "rna"
        }
      }, error = function (err) {
        print(i)
        print(as.character(.tad_gr[i]))
        print(c(did, gid))
        stop(err)
      })

      atac <- setNames(as.numeric(.atac_mat[dar_id == did, -"dar_id"]), names(.atac_mat)[-which(names(.atac_mat) == "dar_id")])

      atac_gene <- atac[match(.valid_samples[[key]], names(atac))]
      gene <- gene[match(.valid_samples[[key]], names(gene))]
      
      if(.randomize) {
        atac_gene <- sample(atac_gene)
        gene <- sample(gene)
      }

      current_cormat[ii,c("pearson","spearman")] <- c(round(cor(atac_gene, gene, method = "pearson"), 2),
                                                      round(cor(atac_gene, gene, method = "spearman"), 2))

      if(gid %in% rownames(.proteins_mat)) {

        protein <- .proteins_mat[gid,]

        atac_protein <- atac[match(.valid_samples[["proteins"]], names(atac))]
        protein <- protein[match(.valid_samples[["proteins"]], names(protein))]
        
        if (.randomize) {
          atac_protein <- sample(atac_protein)
          protein <- sample(protein)
        }

        current_cormat[ii,c("pearson_proteins","spearman_proteins")] <- c(round(cor(atac_protein, protein, method = "pearson"), 2),
                                                                          round(cor(atac_protein, protein, method = "spearman"), 2))
      }

      # TODO: add comparison
      # current_cormat[ii, "comparison"] <- NA

    }
    return(current_cormat)
  }

  dars_in_tad <- findOverlaps(dars_gr, tad_gr)
  dars_in_tad <- cbind(queryHits(dars_in_tad), subjectHits(dars_in_tad))

  genes_in_tad <- findOverlaps(genes_gr, tad_gr)
  genes_in_tad <- cbind(queryHits(genes_in_tad), subjectHits(genes_in_tad))

  valid_samples <- list(rna = intersect(colnames(atac_mat), colnames(gene_mat)),
                        mirna = intersect(colnames(atac_mat), colnames(mirna_mat)),
                        proteins = intersect(colnames(atac_mat), colnames(proteins_mat)))

  if (debug) {
    tad_gr <- tad_gr[sample(length(tad_gr), 2)]
  }

  corlist <- mclapply(seq_along(tad_gr), worker, mc.cores = ncores)

  nulls <- sapply(corlist,is.null)
  while (any(nulls)) {
    nidx <- which(nulls)
    corlist[nidx] <- mclapply(seq_along(tad_gr)[nidx], worker, mc.cores = ncores)
    nulls <- sapply(colist, id.null)
  }

  cormat <- do.call(rbind, corlist)

  if (debug) {
    return(list(corlist=corlist, cormat=cormat))
  } else {
    return(cormat)
  }

}
