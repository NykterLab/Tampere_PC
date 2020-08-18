#' Import path as GenomicRange. Subtract one to starting position.
#'
#' @param path path to bed file.
#' @return gr a GenomicRange object
#' @export
#' @import rtracklayer
#' @examples
#' importGr("some/path/to/file.bed")
#'
importGr <- function (path) {

  gr <- tryCatch({
    rtracklayer::import(path)
  }, error = function(err) {
    gr <- fread(path)
    gr <- makeGRangesFromDataFrame(gr,
                                   keep.extra.columns = TRUE,
                                   starts.in.df.are.0based = TRUE)
    return(gr)
  })

  start(gr) <- start(gr) - 1
  return(gr)
}
