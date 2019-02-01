setwd("/bmt-data/genomics//projects//atac_workdir/")

library(data.table)
library(DESeq2)

ncores <- floor(parallel::detectCores() * 0.5)
setDTthreads(ncores)

getDarTable <- function (dar_table_path) {
  t <- fread(dar_table_path)
  names(t)[1] <- "dar_id_orig"
  t[,dar_id:=gsub(pattern="(chr.*_[0-9]+_[0-9]+)(?:_.*)$", replacement = "\\1", x = dar_id_orig)]
  t
}

getDars <- function (dar_table_path) {
  dar_table <- getDarTable(dar_table_path)
  idx <- sort(dar_rc_mat_norm[dar_table, on=.(dar_id), which=TRUE, nomatch=0])
  rcmat <- dar_rc_mat_norm[idx, ..colsWhitelist]
  rid <- rcmat[,dar_id]
  rcmat[,dar_id:=NULL]
  m <- as.matrix(rcmat)
  rownames(m) <- rid
  counts <- c(
    dars_in_segment=nrow(fread(gsub(x=dar_table_path, pattern="10k_filter.e90ann.bed", replacement="tsv"))),
    dars_10kbpTSS=nrow(dar_table),
    unique_genes=nrow(unique(dar_table[,"Nearest Ensembl"])))
  list(table=dar_table, mat=m, counts=counts)
}


deg_blacklist <- c("BPH_337", "PC_4538", "PC_6102", "PC_6864", "CRPC_531")
dar_ge_corr_blacklist <- c(deg_blacklist, "PC_15194", "PC_22603", "PC_4786")

if(!"rc_rna"%in%ls()){
  load("R_data/RNAdds.RData")

  rna_dds <- rna_dds[,!colnames(rna_dds) %in% dar_ge_corr_blacklist]
  rc_rna  <- assay(rna_dds, "norm.counts")
  
  gr       <- rowRanges(rna_dds)
  gr_sqn   <- as.character(seqnames(gr))
  gr_start <- ifelse(strand(gr) == "+", start(gr), end(gr))

}

if(!"dar_rc_mat_norm"%in%ls()){
  DAR_MAT_PATH <- "dar/RCmat_mor_normed/turp_corrected.tsv.gz"

  dar_rc_mat_norm <- fread(cmd = sprintf("zcat %s", DAR_MAT_PATH), sep ="\t", header = TRUE, showProgress = FALSE)
  names(dar_rc_mat_norm)[1] <- "chrom"
  names(dar_rc_mat_norm) <- gsub(pattern = "^(?:TURP|RP)_(.*)", replacement="\\1", x = names(dar_rc_mat_norm))
  dar_rc_mat_norm[,dar_id := paste(chrom, start, end, sep = "_")]
  cols2keep <- !names(dar_rc_mat_norm) %in% c("chrom", "start", "end")
  dar_rc_mat_norm <- dar_rc_mat_norm[, ..cols2keep]
}

colsWhitelist <- !names(dar_rc_mat_norm) %in% dar_ge_corr_blacklist

columns <- intersect(colnames(rc_rna), names(dar_rc_mat_norm))
columns <- columns[!columns %in% dar_ge_corr_blacklist]

shared.samples <- list(bph  = columns[grepl(pattern = "^BPH", x = columns)],
                       pc   = columns[grepl(pattern = "^PC", x = columns)],
                       crpc = columns[grepl(pattern = "^CRPC", x = columns)])

DAR_FOLDER <- "dar/turp_corrected_log2ratio_2/with_comparison_group_col/homer/"
#dar_files <- list.files(DAR_FOLDER, full.names = TRUE, pattern = "ann.default.bed")
dar_files <- c("dar/turp_corrected_log2ratio_2/with_comparison_group_col/homer//combined_3_comparisons.ann.default.bed")



output.dir <- "dar/correlation/all_genes/"
if(!dir.exists(output.dir)) dir.create(output.dir, recursive=TRUE)

doCor <- function (j, method) {
  rdar <- rc_dar[j,]
  row <- rep(NA,nrow(rc_rna))
  for (i in seq(nrow(rc_rna))) {
    rrna <- rc_rna[i,]
    stopifnot(all(names(rc_rna) == names(rc_dar)))
    
    row[i] <- cor(rrna, rdar, method = method, use = "na.or.complete")
  }
  return(row)
}

for (fp in dar_files) {
  dar_object <- getDars(fp)

  dar_table <- dar_object$table
  rc_dar <- dar_object$mat
  rc_dar <- rc_dar[,match(colnames(rc_rna), colnames(rc_dar))]
  stopifnot(all(colnames(rc_rna) == colnames(rc_dar)))

  
  
  res <- do.call(rbind, mclapply(seq(nrow(rc_dar)), doCor, method = "pearson", mc.cores = ncores))
  dimnames(res) <- list(rownames(rc_dar), rownames(rc_rna))
  write.table(res, file = file.path(output.dir, gsub(x = basename(fp), pattern = "(.*)\\.ann\\.default\\.bed", replacement = "\\1_pearson.tsv")), sep = "\t")
  
  
  
  res <- do.call(rbind, mclapply(seq(nrow(rc_dar)), doCor, method = "spearman", mc.cores = ncores))
  dimnames(res) <- list(rownames(rc_dar), rownames(rc_rna))
  write.table(res, file = file.path(output.dir, gsub(x = basename(fp), pattern = "(.*)\\.ann\\.default\\.bed", replacement = "\\1_spearman.tsv")), sep = "\t")
  
  
  
  res <- do.call(rbind, mclapply(seq(nrow(rc_dar)), function (j) {
    dar_id       <- rownames(rc_dar)[j]
    dar_chr      <- gsub(pattern = "(chr.*)_[0-9]+_[0-9]+", replacement = "\\1", x = dar_id)
    dar_coord    <- as.numeric(c(gsub(pattern =  "chr.*_([0-9]+)_[0-9]+", replacement = "\\1", x = dar_id),
                                 gsub(pattern = "chr.*_[0-9]+_([0-9]+)", replacement = "\\1", x = dar_id)))
    dar_midpoint <- dar_coord[1] + (dar_coord[2] - dar_coord[1]) / 2
    wkeep        <- which(gr_sqn == dar_chr)
    row          <- rep(NA,nrow(rc_rna))
    for (i in seq_along(wkeep)) {
      cur      <- wkeep[i]
      row[cur] <- gr_start[cur] - dar_midpoint
    }
    return(row)
  }, mc.cores = ncores))
  dimnames(res) <- list(rownames(rc_dar), rownames(rc_rna))
  write.table(res, file = file.path(output.dir, gsub(x = basename(fp), pattern = "(.*)\\.ann\\.default\\.bed", replacement = "\\1_distance.tsv")), sep = "\t")

  # stop()

}






