#' computeAtacTSS
#' @description compute log2 fold change, abs med difference, p-value for TSS
#' @param atac_tss Data.table of signals at TSS
#' @export
computeAtacTSS <- function (atac_tss){

  bph_samples <- grep("BPH", names(atac_tss))
  pc_samples <- grep("^PC", names(atac_tss))
  crpc_samples <- grep("^CRPC", names(atac_tss))

  # find medians of each group
  lfc <- data.table(gene_id = atac_tss$gid,
                    bph_med = colMedians(apply(atac_tss[, ..bph_samples, with = FALSE], 1, as.numeric)),
                    pc_med = colMedians(apply(atac_tss[, ..pc_samples, with = FALSE], 1, as.numeric)),
                    crpc_med = colMedians(apply(atac_tss[, ..crpc_samples, with = FALSE], 1, as.numeric)))

  # compute log2 ratio and abs_med_diff
  lfc[, `:=`(
    atac_log2_ratio_pc_bph = round(log2(pc_med) - log2(bph_med), 2),
    atac_log2_ratio_crpc_pc = round(log2(crpc_med) - log2(pc_med), 2),

    atac_abs_med_diff_pc_bph = round(abs(pc_med - bph_med), 2),
    atac_abs_med_diff_crpc_pc = round(abs(crpc_med - pc_med), 2)
  )]

  # Wilcox Mann Withney U-test
  lfc[, `:=`(atac_p_value_pc_bph = NA_real_,
             atac_p_value_crpc_pc = NA_real_)]
  for (i in seq(nrow(lfc))) {
    lfc[i,"atac_p_value_pc_bph"] <- wilcox.test(x = as.numeric(atac_tss[i, ..bph_samples]),
                                                y = as.numeric(atac_tss[i, ..pc_samples]),
                                                exact = FALSE,
                                                paired = FALSE,
                                                alternative = "two.sided")$p.value
    lfc[i,"atac_p_value_crpc_pc"] <- wilcox.test(x = as.numeric(atac_tss[i, ..pc_samples]),
                                                 y = as.numeric(atac_tss[i, ..crpc_samples]),
                                                 exact = FALSE,
                                                 paired = FALSE,
                                                 alternative = "two.sided")$p.value
  }
  lfc[, `:=`(atac_p_value_pc_bph = round(p.adjust(atac_p_value_pc_bph, method = "fdr"), 6),
             atac_p_value_crpc_pc = round(p.adjust(atac_p_value_crpc_pc, method = "fdr"), 6))]
  return(lfc)
}
