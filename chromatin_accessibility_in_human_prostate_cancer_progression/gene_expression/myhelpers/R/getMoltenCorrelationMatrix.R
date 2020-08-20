#' getMoltenCorrelationMatrix
#' @description  Compute Pearson and Spearman correlation matrices between the two input matrices. Melt each and perform an inner join. Return the molten joined matrix.
#' @param m1 first matrix
#' @param m2 second matrix
#' @param varnames character vector. Names to assign to molten columns. First varname should describe m1, second describe m2
#' @param debug set debug mode. Compute correlation only for a small subset of rows of both input matrices.
#' @return data.table
#' @export
#' @examples
#'
getMoltenCorrelationMatrix <- function (m1, m2, varnames, ncores = detectCores(), debug = FALSE, verbose = FALSE) {

  invokeComputeCorrelation <- function (method) {
    if (debug) {
      obj <- computeCorrelation(method, m1, m2, nc = ncores, debug = debug, verbose = verbose)
      res <- obj[["correlation_matrix"]]
    }else {
      res <- computeCorrelation(method, m1, m2, nc = ncores, debug = debug, verbose = verbose)
    }

    if (verbose) {
      cat(sprintf("[getMoltenCorrelationMatrix] class(%s):", method), class(res), "\n")
      cat(sprintf("[getMoltenCorrelationMatrix] dim(%s):", method), dim(res), "\n")
      # if (is.matrix(pearson)) print(pearson[1:10,1:10])
    }
    return(res)
  }

  if (debug) {
    verbose <- TRUE
  }

  if (verbose) {
    cat("[getMoltenCorrelationMatrix] dim(m1):", dim(m1), "\n")
    cat("[getMoltenCorrelationMatrix] dim(m2):", dim(m2), "\n")
  }

  pearson <- invokeComputeCorrelation("pearson")
  spearman <- invokeComputeCorrelation("spearman")

  mpearson <- data.table::melt(pearson, varnames = varnames, value.name = "pearson")
  mspearman <- data.table::melt(spearman, varnames = varnames, value.name = "spearman")

  merged <- merge(mpearson, mspearman, by = varnames, all = TRUE)

  return(merged)
}
