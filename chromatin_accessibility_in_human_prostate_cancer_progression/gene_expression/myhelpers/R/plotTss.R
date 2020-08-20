#' plotTss
#' @description produce plot as in Figure 3C and 3D
#' @param all_cordf
#' @param deg_cordf
#' @param xlim
#' @param ylim
#' @param interesting_genes
#' @param title
#' @import ggplot2
#' @export
plotTss <- function(all_cordf, deg_cordf, xlim = c(-2,2), ylim = c(-10, 10), interesting_genes, title) {

  makeDf <- function (indf) {
    return(data.frame(x = indf$atac_log2FoldChange,
                      y = indf$rna_log2FoldChange,
                      gene_name = indf$gene_name))
  }

  limitDf <- function (indf, xlim, ylim) {
    return(within(indf, {
      x[x < xlim[1]] <- xlim[1]
      x[x > xlim[2]] <- xlim[2]
      y[y < ylim[1]] <- ylim[1]
      y[y > ylim[2]] <- ylim[2]
    }))
  }

  df <- makeDf(all_cordf)
  df_deg <- makeDf(deg_cordf)

  df <- limitDf(df, xlim, ylim)
  df_deg <- limitDf(df_deg, xlim, ylim)

  df$de <- ifelse(all_cordf$gene_name %in% deg_cordf$gene_name, "Yes", "No")
  df$gene_name <- ifelse(df$gene_name %in% interesting_genes & df$de == "Yes", as.character(df$gene_name), "")

  fit1 <- lm( I(y) ~ 0 + x, data = df)
  deg_fit1 <- lm( I(y) ~ 0 + x, data = df_deg)

  p1 <- ggplot(df, aes(x, y, label = gene_name, color = factor(de))) +
    scale_color_manual(values = c("gray60", "red"), name = "Differential expression") +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    background_grid(major = "xy", minor = "none", colour.major = "gray90") +
    geom_point(size = 0.2) +
    geom_point(data = df_deg, mapping = aes(x,y, label = gene_name), size = 0.2, color = "red") +
    geom_text_repel(force = 5, segment.colour = "grey30") +
    geom_abline(slope = coef(fit1), intercept = 0, color = "gray40") +
    geom_abline(slope = coef(deg_fit1), intercept = 0, color = "firebrick") +
    xlab("ATAC-seq log2 fold change") +
    ylab("RNA-seq log2 fold change") +
    ggtitle(title) +
    annotate("text", x = 1.75,  y = 9.5,  size = 5,
             label = sprintf("%d/%d",
                             sum(deg_cordf$atac_log2FoldChange > 0 & deg_cordf$rna_log2FoldChange > 0),
                             sum(all_cordf$quadrant == 1))) +
    annotate("text", x = 1.75,  y = -9.5, size = 5,
             label = sprintf("%d/%d",
                             sum(deg_cordf$atac_log2FoldChange > 0 & deg_cordf$rna_log2FoldChange < 0),
                             sum(all_cordf$quadrant == 2))) +
    annotate("text", x = -1.75, y = -9.5, size = 5,
             label = sprintf("%d/%d",
                             sum(deg_cordf$atac_log2FoldChange < 0 & deg_cordf$rna_log2FoldChange < 0),
                             sum(all_cordf$quadrant == 3))) +
    annotate("text", x = -1.75, y = 9.5,  size = 5,
             label = sprintf("%d/%d",
                             sum(deg_cordf$atac_log2FoldChange < 0 & deg_cordf$rna_log2FoldChange > 0),
                             sum(all_cordf$quadrant == 4))) +
    theme(legend.position = "bottom",
          legend.justification = "center",
          legend.text = element_text(size = 15))

  return(p1)
}
