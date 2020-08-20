#' computeAllDarVsGeneExpressionCorrelations
#' @description perform correlation analysis for each row of atac_mat against each row of gene_mat and mirna_mat. Saves one Rdata file for each atac_mat row in tmp_dir
#' @param genes_gr GRanges object for genes and mirna
#' @param gene_mat Matrix reporting expression values of protein coding genes
#' @param mirna_mat Matrix of expression values of mirna
#' @param atac_mat Matrix of accessibility levels of ATAC-seq regions
#' @param samples_blacklist Blacklisted samples
#' @param tmp_dir Path to tmp dir where store computation results
#' @param mode If "proteins", mirna_mat won't be used
#' @param debug Debug switch
#' @param ncores Number of cores for multiprocessing
#' @return invisible
#' @export
#' @examples
#'
#' @import parallel
computeAllDarVsGeneExpressionCorrelations <- function (genes_gr,
                                                       gene_mat,
                                                       mirna_mat,
                                                       atac_mat,
                                                       samples_blacklist,
                                                       tmp_dir,
                                                       mode = "genes",
                                                       debug = FALSE,
                                                       ncores = detectCores(),
                                                       randomize = FALSE) {

  worker <- function (j, method,
                      .tmp_dir = tmp_dir,
                      .all_gid = all_gid,
                      .atac_mat = atac_mat,
                      .mirna_mat = mirna_mat,
                      .gene_mat = gene_mat,
                      .valid_samples = valid_samples,
                      .mode = mode,
                      .filenames = filenames,
                      .randomize = randomize) {

    dar_id <- rownames(.atac_mat)[j]

    ofp <- file.path(.tmp_dir, method, format(j - (j %% 100), scientific = FALSE), sprintf("%s.Rd", dar_id))

    if (!file.exists(ofp)){

      atac <- .atac_mat[j,]
      row <- setNames(vector("list", length(.all_gid)), .all_gid)

      for (i in seq_along(.all_gid)) {

        gid <- .all_gid[i]

        tryCatch({

          # Try to assign gene expression vector. Two working modes accepted: genes and proteins.
          # If mode == genes and mirna_mat is given, search in mirna_mat or assign from gene_mat.
          # If mode == genes and no mirna_mat is given, search in gene_mat only.
          # If mode == proteins, don't bother searching for mirna, just assign from gene_mat.
          # Raise an error in all other scenarios.
          if (!is.null(mirna_mat) & .mode == "genes"){
            if (gid %in% rownames(.mirna_mat)) {
              gene <- .mirna_mat[gid,]
              key <- "mirna"
            } else {
              gene <- .gene_mat[gid,]
              key <- "rna"
            }
          } else if (.mode == "genes") {
            gene <- .gene_mat[gid,]
            key <- "rna"
          }else if (.mode == "proteins") {
            gene <- .gene_mat[gid,]
            key <- "proteins"
          } else {
            stop("Could not determine gene. Check working mode.")
          }

        }, error = function (err) {

          cat(.mode, gid, gid%in%rownames(.gene_mat), "\n", file = stderr())
          print(err)

        })

        atac_gene <- atac[match(.valid_samples[[key]], names(atac))]
        gene <- gene[match(.valid_samples[[key]], names(gene))]
        
        if(.randomize) {
          atac_gene <- sample(atac_gene)
          gene <- sample(gene)
        }

        row[[i]] <- suppressWarnings(cor.test(atac_gene, gene, method = method, alternative = "two.sided"))
      }

      tryCatch({
        save(row, file=ofp)
      }, error = function (err) {
        print(ofp)
        print(err)
      })
      return(TRUE)
    }
    return(FALSE)
  }

  # doDistance <- function (j,
  #                         .dar_ids = rownames(atac_mat),
  #                         .genes_gr = genes_gr,
  #                         .genes_start = genes_start) {
  #   dar_id <- .dars_ids[j]
  #   dar_chr <- gsub(pattern = "(chr.*)_[0-9]+_[0-9]+", replacement = "\\1", x = dar_id)
  #   dar_coord <- as.numeric(c(start = gsub(pattern =  "chr.*_([0-9]+)_[0-9]+", replacement = "\\1", x = dar_id),
  #                             end = gsub(pattern = "chr.*_[0-9]+_([0-9]+)", replacement = "\\1", x = dar_id)))
  #   dar_midpoint <- (dar_coord[1] + (dar_coord[2] - dar_coord[1])) / 2
  #   wkeep <- which(as.character(seqnames(.genes_gr)) == dar_chr)
  #   row <- rep(NA,length(.genes_gr))
  #   for (i in seq_along(wkeep)) {
  #     cur <- wkeep[i]
  #     row[cur] <- .genes_start[cur] - dar_midpoint
  #   }
  #   return(row)
  # }

  mkdirtree <- function(dir) {
    dir.create(dir, recursive = TRUE)
    for (i in seq(0, nrow(atac_mat))) {
      if (i %% 100 == 0) {
        dir.create(file.path(dir, i))
      }
    }
  }

  ##########
  ## MAIN ##
  ##########

  if (! mode %in% c("proteins", "genes")) {
    stop("Working mode not accepted. Mode must be one of 'genes' or 'proteins'.")
  }

  if (debug & nrow(atac_mat) > 100) {
    unlink(tmp_dir, recursive = TRUE)
    atac_mat <- atac_mat[sample(nrow(atac_mat), 100),]
  }

  if (mode == "genes") {
    valid_samples <- list(rna=intersect(colnames(atac_mat), colnames(gene_mat)),
                          mirna=intersect(colnames(atac_mat), colnames(mirna_mat)))
  } else {
    valid_samples <- list(proteins=intersect(colnames(atac_mat), colnames(gene_mat)))
  }

  if (!is.null(samples_blacklist)) {
    valid_samples <- lapply(valid_samples, function(each) each[!each%in%samples_blacklist])
  }

  if (!is.null(mirna_mat)) {
    all_gid <- c(rownames(gene_mat), rownames(mirna_mat))
  } else {
    all_gid <- rownames(gene_mat)
  }

  if (!dir.exists(tmp_dir)) {
    dir_pearson <- file.path(tmp_dir, "pearson")
    dir_spearman <- file.path(tmp_dir, "spearman")
    mkdirtree(dir_pearson)
    mkdirtree(dir_spearman)
  }

  tryCatch({
    written_files <- simplify2array(mclapply(seq(nrow(atac_mat)), worker, method = "pearson", mc.cores = ncores))
    cat(sprintf("Written %d/%d (%.1f%%) files for pearson computation.\n", sum(written_files), nrow(atac_mat), sum(written_files) / nrow(atac_mat) * 100))
    invisible(gc())
  }, error = function (err) {
    print(head(written_files))
    warnings()
    stop(err)
  })

  tryCatch({
    written_files <- simplify2array(mclapply(seq(nrow(atac_mat)), worker, method = "spearman", mc.cores = ncores))
    cat(sprintf("Written %d/%d (%.1f%%) files for spearman computation.\n", sum(written_files), nrow(atac_mat), sum(written_files) / nrow(atac_mat) * 100))
    invisible(gc())
  }, error = function (err) {
    print(head(written_files))
    warnings()
    stop(err)
  })



  ##########
  # BROKEN #
  ##########

  # genes_start <- ifelse(strand(genes_gr) == "+", start(genes_gr), end(genes_gr))
  # res <- do.call(rbind, mclapply(seq(nrow(atac_mat)), doDistance, mc.cores = ncores))
  # dimnames(res) <- list(rownames(atac_mat), rownames(gene_mat))
  # dist_out <- list(distance_mat=res)
  # save(list = names(dist_out), file = sprintf(fp, "distance"), envir = list2env(dist_out))

  # write.table(res, file = sprintf(fp, "distance"), sep = "\t")

  return(invisible(rownames(atac_mat)))
}
