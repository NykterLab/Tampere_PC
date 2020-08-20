#' computeCorrelation
#' @description  Compute given correlation coefficient matrix between two input matrices. Return matrix of correlations rounded to two decimal places.
#' @param method correlation coefficient calculation method. See `cor` for accepted values
#' @param amat first matrix
#' @param gmat second matrix
#' @param nc number of cores. Defaults to all avaible cores
#' @param debug set debug mode. Compute correlation only for a small subset of rows of both input matrices.
#' @param verbose print log messages
#' @return matrix
#' @import parallel
#' @export
#' @examples
#' m1 <- matrix(rnorm(200), 20, 10, dimnames = list(NULL, letters[1:10]))
#' m2 <- matrix(rnorm(300), 30, 10, dimnames = list(NULL, letters[1:10]))
#' computeCorrelation("pearson", m1,m2)
#'
computeCorrelation <- function (method, amat, gmat,
                                nc = detectCores(),
                                debug = FALSE,
                                verbose = FALSE) {

  if (debug) {
    verbose <- TRUE
    amat <- amat[sample(nrow(amat), 100),]
    gmat <- gmat[sample(nrow(gmat), 100),]
  }

  if (verbose) {
    cat("[computeCorrelation] --INPUTS--\n")
    cat("[computeCorrelation] method is:", method, "\n")
    cat("[computeCorrelation] class(amat):", class(amat), "\n")
    cat("[computeCorrelation] dim(amat):", dim(amat), "\n")
    cat("[computeCorrelation] class(gmat):", class(gmat), "\n")
    cat("[computeCorrelation] dim(gmat):", dim(gmat), "\n")
    cat("[computeCorrelation] -- -- \n")
  }

  intersection <- intersect(colnames(amat), colnames(gmat))

  aidx <- match(intersection, colnames(amat), nomatch = 0)
  gidx <- match(intersection, colnames(gmat), nomatch = 0)

  amat <- amat[,aidx]
  gmat <- gmat[,gidx]

  if (!all(colnames(amat) == colnames(gmat))) {
    cat("colnames(amat):", colnames(amat), "\n")
    cat("colnames(gmat):", colnames(gmat), "\n")
    stop("Error: matrices columns do not match.")
  }

  tryCatch({

    worker <- function(i, .method, .amat, .gmat) {
      atac <- .amat[i,]
      cor_vect <- sapply(setNames(seq(nrow(.gmat)), rownames(.gmat)), function (j, .atac, .meth) {
        gexp <- .gmat[j,]
        cor_value <- round(cor(.atac, gexp, method = .meth), 2)
        return(cor_value)
      },
      .atac = atac,
      .meth = .method)
      return(cor_vect)
    }

    cor_list <- mclapply(setNames(seq(nrow(amat)), rownames(amat)),
                         worker,
                         .method = method,
                         .amat = amat,
                         .gmat = gmat,
                         mc.cores = nc)


    nulls_bool <- sapply(cor_list, is.null)
    while (any(nulls_bool)) {
      nulls <- which(nulls_bool)
      names_nulls <- names(cor_list)[nulls]
      cor_list[names_nulls] <- mclapply(setNames(seq_along(nulls), names_nulls),
                                        worker,
                                        .method = method,
                                        .amat = amat[nulls,],
                                        .gmat = gmat,
                                        mc.cores = nc)
      nulls_bool <- sapply(cor_list, is.null)
    }


    ret <- do.call(rbind, cor_list)

    # dimnames(ret) <- list(rownames(amat), rownames(gmat))

  }, error = function (er) {

    cat("== AN ERROR OCCURRED IN computeCorrelation ==\n")

    cat("--> INPUTS\n")
    cat("method is:", method, "\n")
    cat("class(amat):", class(amat), "\n")
    cat("dim(amat):", dim(amat), "\n")
    cat("class(gmat):", class(gmat), "\n")
    cat("dim(gmat):", dim(gmat), "\n")

    cat("--> OUTPUTS\n")
    cat("class(cor_list):", class(cor_list), "\n")
    cat("length(cor_list):", length(cor_list), "\n")

    nulls <- sapply(cor_list, is.null)
    cat("any(is.null(cor_list)):", any(nulls), "\n")
    if (any(nulls)){
      cat("sum(is.null(cor_list)):", sum(nulls), "\n")
      cat("which(is.null(cor_list)):", names(cor_list)[which(nulls)], "\n")
    }

    if(exists("ret")){
      cat("class(ret):", class(ret), "\n")
      cat("dim(ret):", dim(ret), "\n")
    }

    cat("== ==\n")
    stop(er)
  })

  if (verbose) {
    cat("[computeCorrelation] --OUTPUTS--\n")
    cat("[computeCorrelation] class(cor_list):", class(cor_list), "\n")
    cat("[computeCorrelation] dim(cor_list):", length(cor_list), "\n")
    cat("[computeCorrelation] class(ret):", class(ret), "\n")
    cat("[computeCorrelation] dim(ret):", dim(ret), "\n")

    if (nrow(ret) != nrow(amat)) {
      missing <- rownames(amat)[!rownames(amat) %in% rownames(ret)]
      cat("[computeCorrelation] missing rownames:", missing, "\n")

      for(m in missing) {
        print(cor_list[[m]])
      }

    }

    if (ncol(ret) != nrow(gmat)) {
      missing <- rownames(gmat)[!rownames(gmat) %in% colnames(ret)]
      cat("[computeCorrelation] missing colnames:", missing, "\n")

      for(m in missing) {
        print(cor_list[[m]])
      }

    }

    cat("[computeCorrelation] -- -- \n")
  }

  if (debug) {
    return(list(correlation_matrix=ret, mclapply_output=cor_list))
  } else {
    return(ret)
  }
}

