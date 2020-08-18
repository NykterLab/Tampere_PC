#' dumpFiles
#' @description Write matrix files from files reporting correlations
#' @param targer_dir path to directory holding single Rd files 
#' @param out_fp path to result files (should have 2 placeholders %s)
#' @export
dumpFiles <- function (target_dir, out_fp) {
  
  all_files <- list.files(target_dir, recursive = TRUE, full.names = TRUE)
  
  for (method in c("pearson", "spearman")) {
    
    all_files_matched <- grep(method, all_files, value = TRUE)
    cat(sprintf("%s: detected %d files in %s folder.\n", file.path(target_dir, method), length(all_files_matched), method))
    
    x <- setNames(seq_along(all_files_matched), gsub("(chr.+_[0-9]+_[0-9]+(?:_1)?).+", "\\1", basename(all_files_matched)))
    
    corfile <-  gzfile(sprintf(out_fp, "cormat", method), open = "w")
    pvalfile <-  gzfile(sprintf(out_fp, "pvalue", method), open = "w")
    
    pb <- txtProgressBar(1, length(x), style = 3, width = 40)
    for (i in seq_along(x)) {
      
      fp <- all_files_matched[i]
      dar_id <- names(x)[i]
      
      tryCatch({
        load(fp)
      }, error = function(err) {
        cat(fp,"\n")
        stop(err)
      })
      
      corvect <- setNames(rep(NA, length(row)), names(row))
      pvect <- setNames(rep(NA, length(row)), names(row))
      
      for (j in seq_along(row)) {
        corvect[j] <- round(row[[j]]$estimate, 2)
        pvect[j] <- row[[j]]$p.value
      }
      
      if (i == 1) {
        header <- c("dar_id", names(row))
        cat(header, "\n", file = corfile)
        cat(header, "\n", file = pvalfile)
      }
      
      cat(dar_id, corvect, "\n", file = corfile, append = TRUE)
      cat(dar_id, pvect, "\n", file = pvalfile, append = TRUE)
      
      rm(corvect)
      rm(pvect)
      rm(row)
      invisible(gc())
      
      setTxtProgressBar(pb, i)
    }
    
    close(corfile)
    close(pvalfile)
    
  }
}
