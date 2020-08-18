#' computeDarInTadCorrelationV3
#' @description compute correlation from processed Betools intersect output file
#' @param df the table to process
#' @param mirna_mat matrix of expression values of miRNA
#' @param gene_mat matrix of expression values of genes
#' @param protein_mat matrix of expression values of proteins
#' @param randomize flag. If true shuffle sample order prior of correlation computation
#' @export
computeDarInTadCorrelationV3 <- function (df,
                                          atac_mat,
                                          mirna_mat,
                                          gene_mat,
                                          protein_mat,
                                          samples_blacklist,
                                          id_name,
                                          randomize = FALSE,
                                          ncores = detectCores()) {

  worker <- function(i, .df = df,
                     .atac_mat = atac_mat,
                     .mirna_mat = mirna_mat,
                     .gene_mat = gene_mat,
                     .protein_mat = protein_mat,
                     .samples = samples,
                     .randomize = randomize,
                     .id_name = id_name) {

    gene_id <- .df$gene_id[i]
    peak_id <- unlist(.df[i, ...id_name])

    ret <- setNames(rep(NA_real_, 8), c("pearson", "pearson_pvalue",
                                        "spearman", "spearman_pvalue",
                                        "pearson_protein", "pearson_protein_pvalue",
                                        "spearman_protein", "spearman_protein_pvalue"))


    atac <- .atac_mat[peak_id,]

    if (gene_id %in% rownames(.mirna_mat)) {
      gene <- .mirna_mat[gene_id,]
      key <- "mirna"
    } else {
      gene <- .gene_mat[gene_id,]
      key <- "rna"
    }

    atac <- atac[match(.samples[[key]], names(atac))]
    gene <- gene[match(.samples[[key]], names(gene))]

    if (.randomize) {
      atac <- sample(atac)
      gene <- sample(gene)
    }

    if (sd(atac, na.rm=TRUE) > 0 & sd(gene, na.rm=TRUE) > 0){

      pears <- cor.test(atac,gene)
      spear <- cor.test(atac,gene, method="spearman")

      ret[c("pearson", "pearson_pvalue",
            "spearman", "spearman_pvalue")] <- c(round(pears$estimate, 2), pears$p.value,
                                                 round(spear$estimate, 2), spear$p.value)
    }

    if(gene_id %in% rownames(.protein_mat)) {

      protein <- .protein_mat[gene_id,]
      atac <- .atac_mat[peak_id,]

      if (sd(atac, na.rm=TRUE) > 0 & sd(protein, na.rm=TRUE) > 0){
        protein <- protein[match(samples[["protein"]], names(protein))]
        atac <- atac[match(samples[["protein"]], names(atac))]

        if (.randomize) {
          atac <- sample(atac)
          protein <- sample(protein)
        }

        pears <- cor.test(atac,protein)
        spear <- cor.test(atac,protein, method="spearman")

        ret[c("pearson_protein", "pearson_protein_pvalue",
              "spearman_protein", "spearman_protein_pvalue")] <- c(round(pears$estimate, 2), pears$p.value,
                                                                   round(spear$estimate, 2), spear$p.value)
      }
    }

    return(ret)
  }

  samples <- list(
    mirna = intersect(colnames(atac_mat), colnames(mirna_mat)),
    rna = intersect(colnames(atac_mat), colnames(gene_mat)),
    protein = intersect(colnames(atac_mat), colnames(protein_mat)))

  samples <- lapply(samples, function(x) x[!x %in% samples_blacklist])

  uid <- paste(df$tad_id, unlist(df[,..id_name]), df$gene_id, sep="_")
  x <- setNames(seq(nrow(df)), uid)
  cors <- mclapply(x, worker, mc.cores = ncores)
  cordf <- do.call(rbind, cors)
  cordf <- as.data.table(cordf, keep.rownames = TRUE)

  df$uid <- uid

  ret <- merge(df, cordf, by.x = "uid", by.y="rn")
  ret <- ret[,-1]
  # df <- cbind(df, cormat)
  return(ret)
  # cormat <- do.call(rbind, mclapply(seq(nrow(df)), worker, mc.cores = ncores))
  # df <- cbind(df, cormat)
  # return(df)
}
