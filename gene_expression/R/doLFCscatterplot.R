###########
# HELPERS #
###########

.doLFCscatterCORE <- function (df, ptitle, xlim, ylim, ythr, regression_line, xlab, ylab, show_axis, colors, psize, legend_title, palpha) {

    p <- ggplot(df, aes(x, y)) +
      xlab(xlab) +
      coord_cartesian(xlim = xlim, ylim = ylim) +
      theme_bw() +
      theme(axis.title.y = element_blank(),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14))

    if (show_axis) {
      p <- p +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = 0)
    }

    if (!is.null(ythr)) {
      p <- p + geom_hline(yintercept = c(-ythr, ythr), linetype = "dotted")
    }

    if ("p" %in% colnames(df) & "cor" %in% colnames(df)) {
      p <- p +
          geom_point(aes(size = -log10(p), fill = cor, color = cor), alpha = .4) +
          geom_point(size = psize)
    } else if ("p" %in% colnames(df)) {
      p <- p +
          geom_point(aes(size = -log10(p)), alpha = .4) +
          geom_point(size = psize)
    } else if ("cor" %in% colnames(df)) {
      p <- p +
          geom_point(aes(fill = cor, color = cor), alpha = .4, size = 2) +
          geom_point(size = psize)
    } else if("type" %in% colnames(df)) {
      p <- p + geom_point(aes(color = type), size = psize, alpha = palpha)
    } else {
      p <- p + geom_point(size = psize, alpha = palpha)
    }

    if (regression_line) {
      fit <- lm(y ~ x, data = df)
      r2 <- summary(fit)$r.squared

      xlab <- xlim[2] * 0.85
      ylab <- ylim[2] * 0.95
      lab  <-  sprintf("R2 = %.2f", r2)
      wlab <- strwidth(lab, font = 12, units = "in")
      hlab <- strheight(lab, font = 12, units = "in")
      xlab0 <- xlab - (wlab * 0.5)
      ylab0 <- ylab - (ylab * 0.5)
      #xlab1 <- xlab + (wlab * 0.5)
      #ylab1 <- ylab + (ylab * 0.5)

      if (any(df$x >= xlab0 & df$y >= ylab0, na.rm = T)) {
          # xlab <- xlim[1] * 0.85
        xlab <- xlim[1] * 1.2
      }

      p <- p + geom_smooth(method = "lm") +
          #geom_point(data = data.frame(x = xlab, y = ylab), mapping = aes(x, y), color = "red") +
          annotate("text", x = xlab, y = ylab, label = sprintf("R2 = %.2f", r2))
    }

    if (!is.null(colors) && !is.null(legend_title)) {
      p <- p + scale_color_manual(values=colors, name=legend_title)
    } else if (is.null(colors) && !is.null(legend_title)) {
      p <- p + scale_color_manual(name=legend_title)
    } else if (!is.null(colors) && is.null(legend_title)) {
      p <- p + scale_color_manual(values=colors)
    }


    if ("gene_name" %in% colnames(df)) {
        require("ggrepel")
        p <- p + geom_text_repel(
            data = subset(df, abs(x) >= 0.51 | abs(y) >= 10),
            aes(label = gene_names),
            size = 2,
            box.padding = 0.25,
            point.padding = 0.1,
            segment.color = 'grey70',
            color = 'black')
    }

    return(p)
}

.getAxisLimits <- function (v) {
  c(min(v, na.rm = T), max(v, na.rm = T))
}

.getBinwidth <- function (x) {
  x <- x[!is.na(x)]
  breaks <- pretty(range(x), n = nclass.FD(x), min.n = 1)
  breaks[2] - breaks[1]
}

.doPanel <- function (df, ptitle, xlim, ylim, ythr, regression_line, xlab, ylab, nbins = 100, show_axis, colors, psize, legend_title, palpha, bilayer) {

  show_legend <- "cor" %in% colnames(df) | "p" %in% colnames(df) | "type" %in% colnames(df)


  mainPlot <- .doLFCscatterCORE(df, ptitle, xlim, ylim, ythr, regression_line, xlab, ylab, show_axis, colors, psize, legend_title, palpha)

  if(bilayer) {
    k <- as.numeric(as.factor(df$type)) != 1
    sbst <- df[k,]
    mainPlot <- mainPlot + geom_point(data = sbst, mapping = aes(x,y), color = colors[2],
                                      size = ifelse(length(psize) == 1, psize, psize[k]),
                                      alpha = ifelse(length(palpha) == 1, palpha, palpha[k]))
  }

  if (show_legend) {
    legend <- get_legend(mainPlot)
    mainPlot <- mainPlot + theme(legend.position = "none")
  }

  xbwidth <- .getBinwidth(df$x)
  ybwidth <- .getBinwidth(df$y)

  topHist <- ggplot(df, aes(x)) +
    geom_histogram(binwidth = xbwidth, fill = "grey80", color = "black") +
    scale_x_continuous(limits = xlim) +
    theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14))

  rightHist <- ggplot(df, aes(x = y)) +
    xlab(ylab) +
    geom_histogram(binwidth = ybwidth, fill = "grey80", color = "black") +
    scale_x_continuous(limits = ylim) +
    coord_flip() +
    scale_y_reverse() +
    theme_bw() +
    theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14))

  panel <- plot_grid(NULL, topHist, rightHist, mainPlot,
                     ncol = 2, nrow = 2,
                     align = "hv",
                     rel_widths = c(1, 4), rel_heights = c(1, 1.5))

  if (show_legend) {
    panel <- plot_grid(panel, legend, ncol = 2, rel_widths = c(1, 0.1))
  }

  ptitle_nrow <- length(unlist(strsplit(ptitle, split = "\n"))) + 0.02
  title <- ggdraw() + draw_label(ptitle, fontface = "bold", x = 0.01, hjust = 0)
  panel <- plot_grid(title, panel, nrow = 2, rel_heights = c(ptitle_nrow * 0.08, 1))

  return(panel)
}

########
# MAIN #
########

doLFCscatter <- function (df1,
                          df2 = NULL,
                          ptitle1,
                          ptitle2,
                          main = NULL,
                          outfile = NULL,
                          PICTURES_DIR = "pictures",
                          ythr = 1,
                          xlab = "ATAC LFC",
                          ylab = "RNA LFC",
                          xlim = NULL,
                          ylim = NULL,
                          regression_line = FALSE,
                          nbins = 100,
                          dump_tables = FALSE,
                          tables_path = file.path("tables", "ATACvRNA_TSS"),
                          show_axis = TRUE,
                          return_grob = FALSE,
                          colors = NULL,
                          psize = 0.2,
                          palpha = 1,
                          legend_title = NULL,
                          bilayer = FALSE) {

  stopifnot("x" %in% colnames(df1) & "y" %in% colnames(df1))
  if(!is.null(df2)){ stopifnot("x" %in% colnames(df2) & "y" %in% colnames(df2)) }

  require("ggplot2")
  require("grid")
  require("gridExtra")
  require("cowplot")

  if (is.null(xlim)) {
    xlim <- .getAxisLimits(c(df1$x, df2$x))
  }

  if (is.null(ylim)) {
    ylim <- .getAxisLimits(c(df1$y, df2$y))
  }

  titleGrob <- NULL
  if (!is.null(main)){
    titleGrob <- textGrob(label = main, gp = gpar(fontface = "bold", fontsize = 16))
  }


  topPanel <- .doPanel(df1, ptitle1, xlim, ylim, ythr, regression_line, xlab, ylab, nbins, show_axis, colors, psize, legend_title, palpha, bilayer)

  if(!is.null(df2)) {
    bottomPanel <- .doPanel(df2, ptitle2, xlim, ylim, ythr, regression_line, xlab, ylab, nbins, show_axis, colors, psize, legend_title, palpha, bilayer)
    g <-  arrangeGrob(topPanel, bottomPanel, nrow = 2, top = titleGrob)
  } else {
    g <- arrangeGrob(topPanel, nrow = 1, top = titleGrob)
  }

  if (is.null(outfile)){
    if (return_grob) {
      return(g)
    } else {
      #grid.arrange(topPanel, bottomPanel, nrow = 2, top = titleGrob)
      grid.arrange(g)
    }
  } else {

    if (!dir.exists(PICTURES_DIR)) {
      warnings("Creating ", PICTURES_DIR)
      dir.create(PICTURES_DIR, recursive = T)
    }

    if(!grepl(x = outfile, pattern = "pdf")) {
      outfile <- paste(outfile, "pdf", sep = ".")
    }

    ggsave(g, filename = file.path(PICTURES_DIR, outfile), device = "pdf",
           height = unit(8, "in"), width = unit(10, "in"))

    #ggsave(arrangeGrob(topPanel, bottomPanel, nrow = 2, top = titleGrob),
    #       filename = file.path(PICTURES_DIR, outfile), device = "pdf",
    #       height = unit(8, "in"), width = unit(10, "in"))
  }


  if (dump_tables) {

      #colnames(df1) <- c(xlab, ylab)
      #colnames(df2) <- c(xlab, ylab)

      outfile <- gsub(pattern = "\\.pdf$", "", outfile)

      write.csv(subset(df1, x >= 0 & y >= 0), file = file.path(tables_path, paste(outfile, "TOPsector1.csv", sep = "_")))
      write.csv(subset(df1, x >= 0 & y < 0),  file = file.path(tables_path, paste(outfile, "TOPsector2.csv", sep = "_")))
      write.csv(subset(df1, x < 0 & y < 0),   file = file.path(tables_path, paste(outfile, "TOPsector3.csv", sep = "_")))
      write.csv(subset(df1, x < 0 & y >= 0),  file = file.path(tables_path, paste(outfile, "TOPsector4.csv", sep = "_")))

    if(!is.null(df2)){
      write.csv(subset(df2, x >= 0 & y >= 0), file = file.path(tables_path, paste(outfile, "BOTTOMsector1.csv", sep = "_")))
      write.csv(subset(df2, x >= 0 & y < 0),  file = file.path(tables_path, paste(outfile, "BOTTOMsector2.csv", sep = "_")))
      write.csv(subset(df2, x < 0 & y < 0),   file = file.path(tables_path, paste(outfile, "BOTTOMsector3.csv", sep = "_")))
      write.csv(subset(df2, x < 0 & y >= 0),  file = file.path(tables_path, paste(outfile, "BOTTOMsector4.csv", sep = "_")))
    }

      cat(sprintf("Per-quadrant tables saved in %s.\n", file.path(getwd(), tables_path)))
  }
}

