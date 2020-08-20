quantile.normalize <- function (mat) {
  r <- apply(mat,2,rank,ties.method="first")
  mat_sorted <- apply(mat,2,sort)
  mat_mean <- apply(mat_sorted, 1, mean)
  res <- apply(r,2,function(rnk){ mat_mean[rnk] })
  rownames(res) <- rownames(mat)
  res
}
