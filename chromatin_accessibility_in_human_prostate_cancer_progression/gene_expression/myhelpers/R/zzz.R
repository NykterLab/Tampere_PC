.onLoad <- function (libname, pkgname) {
  suppressPackageStartupMessages(require(GenomicRanges))
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(parallel))
  
  # setDTthreads(detectCores())
}
