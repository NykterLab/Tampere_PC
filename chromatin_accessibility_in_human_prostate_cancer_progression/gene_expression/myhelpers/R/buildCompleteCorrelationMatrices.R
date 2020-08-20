#' buildCompleteCorrelationMatrices
#' @description make complete correlation matrices
#' @param target_dir path to directory containing files to be merged
#' @param out_fp path to RData file to write results to. Has to contain a substitution pattern %s
#' @param proteins if TRUE change names of saved objects including the word "proteins"
#' @param ncores number of cores to use
#' @example
#' @import parallel
#' @export
buildCompleteCorrelationMatrices <- function (target_dir,
                                              out_fp,
                                              proteins = FALSE,
                                              ncores = detectCores()) {

  worker <- function (i, mode, .all_files) {
    tryCatch({
      fp <-.all_files[i]
      dar_id <- gsub("(chr.+_[0-9]+_[0-9]+)_.+", "\\1", basename(fp))
      load(fp)
      rowvect <- setNames(rep(NA, length(row)), names(row))
      for (j in seq_along(row)) {
        cur <- row[[j]]
        if (mode == "cormat") {
          rowvect[j] <- round(cur$estimate, 2)
        } else if (mode == "pvalue") {
          rowvect[j] <- cur$p.value
        }
      }
      return(rowvect)
    }, error = function(err) {
      print(cur)
      print(err)
    })
  }
  
  build <- function (method) {
    
    all_files_matched <- grep(method, all_files, value = TRUE)
    cat(sprintf("%s: detected %d files in %s folder.\n", file.path(target_dir, method), length(all_files_matched), method))
    
    x <- setNames(seq_along(all_files_matched), gsub("(chr.+_[0-9]+_[0-9]+(?:_1)?).+", "\\1", basename(all_files_matched)))
    
    for (mode in c("cormat", "pvalue")) {
      
      fileout <- sprintf(out_fp, mode, method)
      
      outlist <- setNames(vector("list", 1), mode)
      
      if (!file.exists(fileout)) {
        
        outlist[[mode]] <- do.call(rbind,  mclapply(x, worker, mode = mode, .all_files = all_files_matched, mc.cores = ncores))
        
        print(outlist[[mode]][sample(nrow(outlist[[mode]]), 100), sample(ncol(outlist[[mode]]), 100)])
        
        if (proteins) {
          names(outlist) <- paste0(names(outlist), sprintf("_proteins_%s", method)) 
        } else {
          names(outlist) <- paste0(names(outlist), sprintf("_%s", method))
        }
        
        cat(sprintf("Saving %s...", fileout))
        save(list = names(outlist), envir = list2env(outlist), file = fileout)
        cat("\tDONE.\n")
        
        rm(outlist)
        invisible(gc())
        
      } else {
        cat(sprintf("%s exists.\n",fileout))
      }
    }
    
    # full_cormat <- do.call(rbind,  mclapply(x, worker, mode = "cor", .all_files = all_files_matched, mc.cores = ncores))
    # full_pvmat <- do.call(rbind, mclapply(x, worker, mode = "pvalue", .all_files = all_files_matched, mc.cores = ncores))
    # 
    # outlist <- list(full_cormat, full_pvmat)
    # 
    # if (proteins) {
    #   names(outlist) <- c(sprintf("cormat_proteins_%s", method), sprintf("pvalue_proteins_%s", method))
    # } else {
    #   names(outlist) <- c(sprintf("cormat_%s", method), sprintf("pvalue_%s", method))
    # }
    # 
    # save(list = names(outlist), envir = list2env(outlist), file = sprintf(out_fp, method))
    # 
    # rm(full_cormat, full_pvmat)
    # invisible(gc())
  }

  ##########
  ## MAIN ##
  ##########

  all_files <- file.path(target_dir, list.files(target_dir, recursive = TRUE))

  build("pearson")
  build("spearman")

  invisible()
}
