#' mergeGtrd
#' @description Find TFBS overalpping features
#' @param gr ranges to search TFBS for
#' @param gtrd GTRD GenomicRanges
#' @return data.table
#' @export
#' @examples
#' mergeGtrd(dt, gtrd)
#' @import GenomicRanges
#' @import data.table
mergeGtrd <- function (gr, gtrd) {

  # df <- data.frame(chr   = gsub(pattern = "(chr.*)_[0-9]+_[0-9]+", replacement = "\\1", x = dt[["dar_id"]]),
  #                  start = gsub(pattern = "chr.*_([0-9]+)_[0-9]+", replacement = "\\1", x = dt[["dar_id"]]),
  #                  end   = gsub(pattern = "chr.*_[0-9]+_([0-9]+)", replacement = "\\1", x = dt[["dar_id"]]),
  #                  name  = dt[["dar_id"]],
  #                  stringsAsFactors = FALSE)

  # gr    <- makeGRangesFromDataFrame(df, starts.in.df.are.0based = FALSE, keep.extra.columns = TRUE, ignore.strand = TRUE)

    olaps <- findOverlaps(gtrd, gr)

  tf_dt <- data.table(name = names(gr)[subjectHits(olaps)], gtrd_tf = gtrd$name[queryHits(olaps)])
  tf_dt <- tf_dt[, .(gtrd_tf = paste(unique(gtrd_tf), collapse = ", ")), by = name]
  setnames(tf_dt, c("name", "gtrd"))

  return(tf_dt)
}
