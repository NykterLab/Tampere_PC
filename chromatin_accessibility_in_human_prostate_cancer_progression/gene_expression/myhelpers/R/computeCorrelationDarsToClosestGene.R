#' computeCorrelationDarToClosestGene
#' @description compute correlation dar to closest gene according to Homer annotation
#' @param dar_to_closest_gene data.table output from Homer
#' @param mirna_map mirBase to Ensembl gene id mapping
#' @param mirna_mat matrix of miRNA expression values
#' @param gene_mat matrix of genes expression values
#' @param atac_mat matrix of atac accessibility values
#' @param proteins_mat matrix of protein expression values
#' @param randomize flag. If true shuffle samples order prior to correlation computation
#' @export
computeCorrelationDarToClosestGene <- function (dar_to_closest_gene, 
                                                mirna_map, 
                                                mirna_mat, 
                                                gene_mat, 
                                                atac_mat, 
                                                protein_mat, 
                                                samples_blacklist,
                                                id_column,
                                                randomize = FALSE,
                                                ncores = detectCores()) {
  
  worker <- function (i, 
                      .dar_to_closest_gene = dar_to_closest_gene, 
                      .mirna_map = mirna_map, 
                      .mirna_mat = mirna_mat, 
                      .gene_mat = gene_mat,
                      .atac_mat = atac_mat,
                      .protein_mat = protein_mat,
                      .samples = samples,
                      .randomize = randomize,
                      .id_column = id_column,
                      .names_vector = names_vector) {
    
    dar_id <- unlist(.dar_to_closest_gene[i, ...id_column])
    gene_id <- .dar_to_closest_gene$`Nearest Ensembl`[i]
    ret <- setNames(rep(NA,8), .names_vector)
    
    if (gene_id %in% .mirna_map$gene_id) {
      gene_id <- .mirna_map[.mirna_map$mirna_id, "mirna_name"]
      gene <- .mirna_mat[gene_id,]
      key <- "mirna"
    } else if (gene_id %in% rownames(.gene_mat)) {
      gene <- .gene_mat[gene_id,]
      key <- "rna"
    } else {
      return(ret)
    }
    
    tryCatch({
      atac <- .atac_mat[dar_id,]  
    }, error = function (er) {
      print(dar_id)
      return(er)
    })
    
    
    atac <- atac[match(.samples[[key]], names(atac))]
    gene <- gene[match(.samples[[key]], names(gene))]
    
    if (.randomize) {
      atac <- sample(atac)
      gene <- sample(gene)
    }
    
    if (sd(atac, na.rm=TRUE) > 0 & sd(gene, na.rm=TRUE) > 0) {
      pears <- cor.test(atac, gene)
      spear <- cor.test(atac, gene, method = "spearman")
      ret[c("pearson", "pearson_pvalue", 
            "spearman", "spearman_pvalue")] <- c(round(pears$estimate,2), pears$p.value, 
                                                 round(spear$estimate,2), spear$p.value)
    }
    
    if(gene_id %in% rownames(.protein_mat)) {
      
      protein <- .protein_mat[gene_id,]
      atac <- .atac_mat[dar_id,]
      
      protein <- protein[match(.samples[["proteins"]], names(protein))]
      atac <- atac[match(.samples[["proteins"]], names(atac))]
      
      if (.randomize) {
        protein <- sample(protein)
        atac <- sample(atac)
      }
      
      if (sd(atac, na.rm=TRUE) > 0 & sd(protein, na.rm=TRUE) > 0){
        pears <- cor.test(atac, protein)
        spear <- cor.test(atac, protein, method="spearman")
        ret[c("pearson_protein", "pearson_protein_pvalue", 
              "spearman_protein", "spearman_protein_pvalue")] <- c(round(pears$estimate,2), pears$p.value, 
                                                                   round(spear$estimate,2), spear$p.value)
      }
    }
    return(ret)
  }
  
  names_vector <- c("pearson", "pearson_pvalue", 
                    "spearman", "spearman_pvalue", 
                    "pearson_protein", "pearson_protein_pvalue", 
                    "spearman_protein", "spearman_protein_pvalue")
  
  samples <- list(mirna = intersect(colnames(atac_mat), colnames(mirna_mat)), 
                  rna = intersect(colnames(atac_mat), colnames(gene_mat)),
                  proteins = intersect(colnames(atac_mat), colnames(protein_mat)))
  
  samples <- lapply(samples, function(x) x[!x %in% samples_blacklist])
  
  cormat <- t(as.data.table(mclapply(seq(nrow(dar_to_closest_gene)), worker, mc.cores = ncores)))
  cormat <- as.data.table(cormat)
  setnames(cormat, names(cormat), names_vector)
  
  dar_to_closest_gene <- cbind(dar_to_closest_gene, cormat)
  return(dar_to_closest_gene)
  
}
