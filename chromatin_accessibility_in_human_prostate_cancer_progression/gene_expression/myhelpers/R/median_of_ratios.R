# m is a matrix with one row per gene and one column per sample.
mor <- function (m, threshold = 50, 
                 add_one = TRUE,
                 return_coefficients = FALSE, 
                 use_median = FALSE,
                 do_rounding = TRUE,
                 digits = 5
                ) {
  
  if(any(m<0, na.rm = TRUE)) stop("Negative counts detected. Aborting.", call. = FALSE)

  ## shift all cells by 1.
  if(add_one) {
    cat("Adding 1 to all cells. Turn off with add_one=FALSE.\n", file = stderr())
    m <- m + 1
  }
  
  ## build consensus sample
  if(use_median){
    consensus <- apply(m, 1, median, na.rm = TRUE)
  }else{
    consensus <- exp(rowMeans(log(m), na.rm = TRUE)) # geometric mean          
  }

  ## compute median-of-ratios  
  keep <- consensus > threshold
  mors <- apply(m, 2, function (sample) { median(sample[keep] / consensus[keep], na.rm = TRUE) })
    
  ## normalize input
  norm <- sapply(1:ncol(m), function (i) { m[,i] / mors[i] })
  
  ## round
  if (do_rounding) {
    norm <- round(norm, digits = digits)
    if(return_coefficients) mors <- round(mors, digits = digits)
  }
  
  ## keep colnames and rownames
  colnames(norm) <- colnames(m)
  rownames(norm) <- rownames(m)
  
  if(return_coefficients){
    return(list(medians = consensus, 
                coeff = mors, 
                norm = norm))
  }else{
    return(norm)    
  } 
}
