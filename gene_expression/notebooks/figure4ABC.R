#load("/bmt-data/genomics/projects/atac_workdir/gene_expression_correlation/v7/env.rdata")

load(file = "/bmt-data/genomics/projects/atac_workdir/gene_expression_correlation/tf_analysis/counts_and_conmat.rdata")

library(pheatmap)
library(VennDiagram)

# FIGURE 4A

do_panel_A <- function (mat, path,
                        clustering_method = "average",
                        clustering_distance = "manhattan",
                        thr = 100){

    k1 <- apply(mat, 1, function(x) any(x >= thr))
    message(sum(k1))
    mat_filter <- mat[k1, k1]

    pheatmap(mat_filter,
             clustering_method=clustering_method,
             clustering_distance_rows=clustering_distance,
             clustering_distance_cols=clustering_distance,
             show_rownames=FALSE, show_colnames=FALSE,
             main=sprintf("%d/%d TFs - at least one connection >=%d per row\ndistance = %s - clustering = %s", sum(k1), nrow(mat), thr, clustering_distance, clustering_method),
             filename=path,
             width=8, height=8,
             border_color=NA)
    return(mat_filter)

}




# FIGURE 4B

do_panel_B <- function (mat, path,
                        clustering_method = "complete",
                        clustering_distance = "euclidean",
                        thr = 300,
                        hc = NULL
                        ) {
    print(dim(mat))

    k2 <- apply(mat, 1, function(x) any(x >= thr))
    print(sum(k2))
    mat_filter <- mat[k2, k2]
    print(dim(mat_filter))

    if (!is.null(hc)) {
        message("Using hierarchical clustering from input paramters. Ignoring clustering_method and clustering_distance.")
        pheat <- pheatmap(mat_filter, border_color=NA,
              cluster_cols=hc,
              cluster_rows=hc,
              main=sprintf("%d TFs - at least one connection >=%d per row\ndistance = %s - clustering = %s", sum(k2), thr, clustering_distance, clustering_method),
              filename=path,
              width=8, height=8)
    } else {
        pheat <- pheatmap(mat_filter, border_color=NA,
              clustering_distance_rows=clustering_distance,
              clustering_distance_cols=clustering_distance,
              clustering_method=clustering_method,
              main=sprintf("%d TFs - at least one connection >=%d per row\ndistance = %s - clustering = %s", sum(k2), thr, clustering_distance, clustering_method),
              filename=path,
              width=8, height=8)
    }


    return(pheat)

}



# FIGURE 4C
get_genes <- function(cl, target_group, frac = 0.75){


    genes_cl <- unlist(lapply(cl, function(x) genelist[[x]][[target_group]]))
    genes_cl <- unique(genes_cl)

    m1 <- matrix(0, length(genes_cl), length(cl), dimnames=list(genes_cl, cl))
    for (tf in cl) {
        g <- genelist[[tf]][[target_group]]
        m1[g, tf] <- 1
    }

    if (frac > 0){
        rs <- rowSums(m1)
        thr <- ceiling(ncol(m1) * frac)
        message("threshold: ", thr, " / ", ncol(m1))
        k <- rs >= thr
        message("number of genes over threshold: ", sum(k), " / ", length(unique(genes_cl)))
        ret <- rownames(m1)[k]
    } else {
        message("Using all genes")
        ret <- rownames(m1)
    }

    return(ret)

}

do_panel_C <- function(pheat, path, thr = 0.2, target_group = "dars") {
    message(sprintf("Using genes regulated by %d %% of TFs in clusters", thr * 100))

    tr <- pheat$tree_row
    ct <- cutree(tr, 2)

    cl1 <- names(ct[ct == 1])
    cl2 <- names(ct[ct == 2])


    message(sprintf("cluster 1 (%s): %d - cluster 2 (%s): %d",
                    ifelse("AR" %in% cl1, "AR", "other"),
                    length(cl1),
                    ifelse("AR" %in% cl2, "AR", "other"),
                    length(cl2)))


    genes_cl1 <- get_genes(cl1, target_group = target_group, frac = thr)
    genes_cl2 <- get_genes(cl2, target_group = target_group, frac = thr)

    venn.diagram(x=list(AR_cluster=genes_cl1, other=genes_cl2),
                 height=8, width=8, units="in",
                 main=sprintf("AR cluster = %d TFs - other cluster = %d TFs\nGenes regulated by at least %d %% of TFs in cluster.",
                              ifelse("AR" %in% cl1, length(cl1), length(cl2)), # AR cluster
                              ifelse("AR" %in% cl1, length(cl2), length(cl1)), # other cluster
                              thr * 100),
                 cex=3,
                 filename=path, imagetype="svg",
                 fill=c("red", "green"), alpha=0.6,
                 euler.d=TRUE, scaled=TRUE)

}

panel_B_clustering_method = "average"

conmat_filter <- do_panel_A(conmat, "/home/ft413468/test/figure4A.pdf")
pheat <- do_panel_B(conmat_filter, "/home/ft413468/test/figure4B.pdf", clustering_method = panel_B_clustering_method)
do_panel_C(pheat, "~/test/figure4C.svg")

### TEST WITH PEAKS
conmat_peaks_filter <- do_panel_A(conmat_peaks, "/home/ft413468/test/figure4A_peaks.pdf", thr = 1000)
pheat_peaks <- do_panel_B(conmat_peaks_filter, "/home/ft413468/test/figure4B_peaks.pdf", thr = 3000, clustering_method = panel_B_clustering_method)
do_panel_C(pheat_peaks, path = "~/test/figure4C_peaks.svg", target_group = "peaks")

### GET DARs CLUSTERING AND PEAKS GENESET
#do_panel_B(conmat_peaks_filter, "/home/ft413468/test/figure4B_peaks_dar_clustering_upgma.pdf", thr = 3000, hc = pheat$tree_row)
do_panel_C(pheat, path = "~/test/figure4C_peaks_dars_clusters_upgma.svg", target_group = "peaks")
