#' setThreads
#' @description set number of threads used by data.table
#' @param n number of threads
#' @export
setThreads <- function (n) {
  setDTthreads(n)
}
