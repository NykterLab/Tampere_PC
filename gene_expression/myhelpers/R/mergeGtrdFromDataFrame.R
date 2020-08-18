#' MergeGtrdFromDataFrame
#' @description
#' @param input_df data.frame to be merged
#' @param gid_df data frame of coordinates to look for TFBS
#' @param gtrd GenomicRanges object of TFBS coordinates
#' @return data.table
#' @export
#' @examples
#'
mergeGtrdFromDataFrame <- function (input_df, gid_df, gtrd, by) {

  if (class(gid_df)[1] == "GRanges") {
    gr <- gid_df
  } else {
    gr <- makeGRangesFromDataFrame(as.data.frame(gid_df),
                                   starts.in.df.are.0based = FALSE,
                                   keep.extra.columns = TRUE)
    names(gr) <- gr$gid
  }

  gtrd_colname <- gsub("_", " ", deparse(substitute(gtrd)))
  tfbs <- myhelpers::mergeGtrd(gr, gtrd)
  setnames(tfbs, c("name", gtrd_colname))

  ret <- merge(input_df, tfbs, by.x = by, by.y = "name", all.x = TRUE)

  return(ret)
}
