load("/bmt-data/genomics/projects/atac_workdir/gene_expression_correlation/v7/env.rdata")

t <- data.frame(
  tf = unique(gtrd_all$name),
  bs_dars_pc_bph_opening = 0,
  bs_dars_pc_bph_closing = 0,
  bs_dars_crpc_pc_opening = 0,
  bs_dars_crpc_pc_closing = 0,
  bs_peaks = 0,
  dars_corr_genes = 0,
  peaks_corr_genes = 0
)

bs_tot <- read.table("/home/ft413468/projects/atacseq_tamperepc/gtrd_counts")
colnames(bs_tot) <- c("tot_bs_gtrd", "tf")
t <- merge(t, bs_tot, by = "tf")

makeTable <- function (comparison, direction) {
  ret <- rbind(subset(dar_in_prom, subset = status == direction & comparison_group == comparison & (abs(pearson) > 0.5 | abs(spearman) > 0.5),
                      select = c("gene_id", "dar_id", "gtrd all")),
               subset(dar_in_tad, subset = status == direction & comparison_group == comparison & (abs(pearson) > 0.5 | abs(spearman) > 0.5),
                      select = c("gene_id", "dar_id", "gtrd all")))
  x <- subset(dar_to_closest_gene, subset = status == direction & comparison_group == comparison & (abs(pearson) > 0.5 | abs(spearman) > 0.5),
              select = c("map_id", "dar_id", "gtrd all"))
  colnames(x) <- c("gene_id", "dar_id", "gtrd all")
  ret <- rbind(ret, x)
  ret <- ret[!duplicated(ret),]
  ret <- ret[, .(genes = paste(gene_id, collapse = ", "), ngenes = .N), by = .(dar_id, `gtrd all`)]
  return(ret)
}

makePeaksTable <- function() {
  ret <- rbind(subset(peaks_in_prom, abs(pearson) > 0.5 | abs(spearman) > 0.5 , select = c("peak_id", "gene_id", "gtrd all")),
               subset(peaks_in_tad, abs(pearson) > 0.5 | abs(spearman) > 0.5 , select = c("peak_id", "gene_id", "gtrd all")))
  x <- subset(peaks_to_closest_gene, abs(pearson) > 0.5 | abs(spearman) > 0.5 , select = c("peak_id", "map_id", "gtrd all"))
  colnames(x) <- c("peak_id", "gene_id", "gtrd all")
  ret <- rbind(ret, x)
  ret <- ret[!duplicated(ret),]
  ret <- ret[, .(genes = paste(gene_id, collapse = ", "), ngenes = .N), by = .(peak_id, `gtrd all`)]
  return(ret)
}

#makePeaksTable2 <- function() {
#  ret <- rbind(subset(peaks_in_prom, select = c("peak_id", "gene_id", "gtrd all")),
#               subset(peaks_in_tad, select = c("peak_id", "gene_id", "gtrd all")))
#  x <- subset(peaks_to_closest_gene, select = c("peak_id", "map_id", "gtrd all"))
#  colnames(x) <- c("peak_id", "gene_id", "gtrd all")
#  ret <- rbind(ret, x)
#  ret <- ret[!duplicated(ret),]
#  ret <- ret[, .(genes = paste(gene_id, collapse = ", "), ngenes = .N), by = .(peak_id, `gtrd all`)]
#  return(ret)
#}

# all_peaks2 <- makePeaksTable2()
# tf_peaks2 <- strsplit(all_peaks2$`gtrd all`, ", ")
# v <- sapply(tf_peaks2, length)
# pdf("/bmt-data/genomics/projects/atac_workdir/gene_expression_correlation/tf_analysis/tfbs_per_peak.pdf")
# par(mfrow=c(1,2), pty="s")
# hist(v,"fd", main = "Number of TBFS per peak (all peaks)", xlab = "Number of binding sites")

all_peaks <- makePeaksTable()
tf_peaks <- strsplit(all_peaks$`gtrd all`, ", ")
# v <- sapply(tf_peaks, length)
# hist(v,"fd", main = "Number of TBFS per peak (correlating peaks)", xlab = "Number of binding sites")
# dev.off()

genelist <- setNames(vector("list", nrow(t)), t$tf)
genelist <- lapply(genelist, function(each) list(dars=NULL, peaks=NULL, pc_bph=NULL, crpc_pc=NULL))

dars_pc_bph_open <- makeTable("pc_bph", "opening")
tf_pc_bph_open <- strsplit(dars_pc_bph_open$`gtrd all`, ", ")

dars_pc_bph_close <- makeTable("pc_bph", "closing")
tf_pc_bph_close <- strsplit(dars_pc_bph_close$`gtrd all`, ", ")

dars_crpc_pc_open <- makeTable("crpc_pc", "opening")
tf_crpc_pc_open <- strsplit(dars_crpc_pc_open$`gtrd all`, ", ")

dars_crpc_pc_close <- makeTable("crpc_pc", "closing")
tf_crpc_pc_close <- strsplit(dars_crpc_pc_close$`gtrd all`, ", ")


my_grep <- function (each, pat) any(grepl(pat, each))

for (i in seq(nrow(t))) {
# for (i in seq(10)) {
  tf <- t$tf[i]
  pat <- sprintf("^%s$", tf)

  pc_bph_open_with_bs <- sapply(tf_pc_bph_open, my_grep, pat = pat)
  pc_bph_close_with_bs <- sapply(tf_pc_bph_close, my_grep, pat = pat)
  crpc_pc_open_with_bs <- sapply(tf_crpc_pc_open, my_grep, pat = pat)
  crpc_pc_close_with_bs <- sapply(tf_crpc_pc_close, my_grep, pat = pat)

  peaks_with_bs <- sapply(tf_peaks, my_grep, pat = pat)

  t[i, "bs_dars_pc_bph_opening"] <- sum(pc_bph_open_with_bs)
  t[i, "bs_dars_pc_bph_closing"] <- sum(pc_bph_close_with_bs)
  t[i, "bs_dars_crpc_pc_opening"] <- sum(crpc_pc_open_with_bs)
  t[i, "bs_dars_crpc_pc_closing"] <- sum(crpc_pc_close_with_bs)
  t[i, "bs_dars"] <- t[i, "bs_dars_pc_bph_opening"] + t[i, "bs_dars_pc_bph_closing"] + t[i, "bs_dars_crpc_pc_opening"] + t[i, "bs_dars_crpc_pc_closing"]

  t[i, "bs_peaks"] <- sum(peaks_with_bs)

  genes_pc_bph_open <- strsplit(dars_pc_bph_open$genes[pc_bph_open_with_bs], ", ")
  genes_pc_bph_close <- strsplit(dars_pc_bph_close$genes[pc_bph_close_with_bs], ", ")
  genes_crpc_pc_open <- strsplit(dars_crpc_pc_open$genes[crpc_pc_open_with_bs], ", ")
  genes_crpc_pc_close <- strsplit(dars_crpc_pc_open$genes[crpc_pc_close_with_bs], ", ")

  genes_dars <- Reduce(union, c(genes_pc_bph_open, genes_pc_bph_close, genes_crpc_pc_open, genes_crpc_pc_close))
  t[i, "dars_corr_genes"] <- length(genes_dars)

  genes_pc_bph <- Reduce(union, c(genes_pc_bph_open, genes_pc_bph_close))
  genes_crpc_pc <- Reduce(union, c(genes_crpc_pc_open, genes_crpc_pc_close))

  genes_peaks <- strsplit(all_peaks$genes[peaks_with_bs], ", ")
  genes_peaks <- Reduce(union, genes_peaks)
  t[i, "peaks_corr_genes"] <- length(genes_peaks)

  genelist[[tf]][["dars"]] <- genes_dars
  genelist[[tf]][["peaks"]] <- genes_peaks

  genelist[[tf]][["pc_bph"]] <- genes_pc_bph
  genelist[[tf]][["crpc_pc"]] <- genes_crpc_pc

}

conmat <- matrix(0, nrow(t), nrow(t), dimnames = list(t$tf, t$tf))
conmat_pc_bph <- matrix(0, nrow(t), nrow(t), dimnames = list(t$tf, t$tf))
conmat_crpc_pc <- matrix(0, nrow(t), nrow(t), dimnames = list(t$tf, t$tf))
conmat_peaks <- matrix(0, nrow(t), nrow(t), dimnames = list(t$tf, t$tf))

for (i in seq(nrow(conmat)-1)) {
  tfi <- rownames(conmat)[i]

  gi <- genelist[[tfi]][["dars"]]
  gi_pc_bph <- genelist[[tfi]][["pc_bph"]]
  gi_crpc_pc <- genelist[[tfi]][["crpc_pc"]]
  gi_peaks <-genelist[[tfi]][["peaks"]]

  for (j in seq(i+1, ncol(conmat))) {
    tfj <- colnames(conmat)[j]

    gj <- genelist[[tfj]][["dars"]]
    gj_pc_bph <- genelist[[tfj]][["pc_bph"]]
    gj_crpc_pc <- genelist[[tfj]][["crpc_pc"]]
    gj_peaks <-genelist[[tfj]][["peaks"]]

    conmat[i,j] <- length(intersect(gi,gj))
    conmat[j,i] <- conmat[i,j]

    conmat_pc_bph[i,j] <- length(intersect(gi_pc_bph, gj_pc_bph))
    conmat_pc_bph[j,i] <- conmat_pc_bph[i,j]

    conmat_crpc_pc[i,j] <- length(intersect(gi_crpc_pc, gj_crpc_pc))
    conmat_crpc_pc[j,i] <- conmat_crpc_pc[i,j]

    conmat_peaks[i,j] <- length(intersect(gi_peaks, gj_peaks))
    conmat_peaks[j,i] <- conmat_peaks[i,j]

  }
}
diag(conmat) <- 0
diag(conmat_pc_bph) <- 0
diag(conmat_crpc_pc) <- 0
diag(conmat_peaks) <- 0

save(t, conmat, conmat_pc_bph, conmat_crpc_pc, conmat_peaks, genelist,
     file = "/bmt-data/genomics/projects/atac_workdir/gene_expression_correlation/tf_analysis/counts_and_conmat.rdata")

load(file = "/bmt-data/genomics/projects/atac_workdir/gene_expression_correlation/tf_analysis/counts_and_conmat.rdata")


tf <- rownames(conmat)
tf <- merge(as.data.frame(tf), ann.gene, by.x="tf", by.y="gene_name", all.x=TRUE)
tf <- subset(tf, select=c("tf", "gene_id"))
# rna_tf <- rna_dds[tf$gene_id,]
rna_tf <- rna_dds[rownames(rna_dds) %in% tf$gene_id,]
unamapped <- tf$gene_id[which(!tf$gene_id %in% rownames(rna_dds))]
#           tf         gene_id
# 65     CSHL1 ENSG00000204414
# 122   FERD3L ENSG00000146618
# 152    FOXR2 ENSG00000189299
# 495    SOX10 ENSG00000100146
# 585    ZBED1 ENSG00000214717
# 686  ZNF280A ENSG00000169548
# 842 ZSCAN5DP ENSG00000267908

med <- apply(assay(rna_tf, "norm.counts"), 1, median)
# q25 <- quantile(med,0.25)
t <- 4.22
expressed_ensg <- names(med)[med > t]
expressed_tf <- tf[tf[,2] %in% expressed_ensg, 1]

conmat <- conmat[rownames(conmat) %in% expressed_tf, colnames(conmat) %in% expressed_tf]
conmat_peaks <- conmat_peaks[rownames(conmat_peaks) %in% expressed_tf, colnames(conmat_peaks) %in% expressed_tf]
conmat_pc_bph <- conmat_pc_bph[rownames(conmat_pc_bph) %in% expressed_tf, colnames(conmat_pc_bph) %in% expressed_tf]
conmat_crpc_pc <- conmat_crpc_pc[rownames(conmat_crpc_pc) %in% expressed_tf, colnames(conmat_crpc_pc) %in% expressed_tf]


library(pheatmap)

generate_heatmap <- function(mat, level, main, root_path, branch = NULL) {
   message(level, branch, ": ", nrow(mat), " rows - ", root_path)
   if (!is.matrix(mat)) {
       return(NULL)
   }
   if (level == 7 | nrow(mat) < 5) {
        return(NULL)
   }

   dir.create(root_path, recursive=TRUE, showWarnings=FALSE)

   filename <- ifelse(is.null(branch), sprintf("conmat%d.png", level), sprintf("conmat%d%s.png", level, branch))
   fp <- file.path(root_path,filename)

   tablename <- ifelse(is.null(branch), sprintf("conmat%d.csv", level), sprintf("conmat%d%s.csv", level, branch))
   write.csv(subset(t, subset = tf %in% rownames(mat)), file=file.path(root_path, tablename))

   rdsname <- ifelse(is.null(branch), sprintf("conmat%d.rds", level), sprintf("conmat%d%s.rds", level, branch))
   saveRDS(mat, file = file.path(root_path, rdsname))

   ph <- pheatmap(mat,
                  main = sprintf("%s - %d TF", main, nrow(mat)),
                  clustering_method = "complete",
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  show_rownames = nrow(mat) < 50,
                  show_colnames = nrow(mat) < 50,
                  filename = file.path(root_path,filename))

   ct <- cutree(ph$tree_row, 2)

   cluster1 <- names(ct)[ct == 1]
   mat1 <- mat[cluster1, cluster1]
   message("Cluster 1: ", length(cluster1))

   cluster2 <- names(ct)[ct == 2]
   mat2 <- mat[cluster2, cluster2]
   message("Cluster 2: ", length(cluster2))

   next_level <- level + 1

   if (length(cluster1) > 0) {
       branch <- "A"
       generate_heatmap(mat1, next_level, main, root_path = file.path(root_path, paste0(level, branch)), branch = branch)
   }

   if (length(cluster2) > 0) {
       branch <- "B"
       generate_heatmap(mat2, next_level, main, root_path = file.path(root_path, paste0(level, branch)), branch = branch)
   }

   return(NULL)
}

generate_heatmap(conmat, 0, "ALL DARs", "/bmt-data/genomics/projects/atac_workdir/gene_expression_correlation/tf_analysis_expressed_tf_only/all_dars")
generate_heatmap(conmat_pc_bph, 0, "PCvsBPH", "/bmt-data/genomics/projects/atac_workdir/gene_expression_correlation/tf_analysis_expressed_tf_only/pc_bph")
generate_heatmap(conmat_crpc_pc, 0, "CRPCvsPC", "/bmt-data/genomics/projects/atac_workdir/gene_expression_correlation/tf_analysis_expressed_tf_only/crpc_pc")
generate_heatmap(conmat_peaks, 0, "Peaks", "/bmt-data/genomics/projects/atac_workdir/gene_expression_correlation/tf_analysis_expressed_tf_only/peaks")
#generate_heatmap(conmat, 0, "ALL DARs", "/bmt-data/genomics/projects/atac_workdir/gene_expression_correlation/tf_analysis/all_dars")


# generate_heatmap2 <- function(mat, level) {
#     message(level)
#     if (!is.matrix(mat)) {
#        return(NULL)
#     }
#     if(level == 10 | nrow(mat) < 5){
#         return(NULL)
#     }
#     d <- dist(mat, method="euclidean")
#     hc <- hclust(d, method = "complete")
#     ct <- cutree(hc, 2)
#     cl1 <- names(ct)[ct==1]
#     cl2 <- names(ct)[ct==2]
#     return(list(level=level,
#                 members=rownames(mat),
#                 branch1=generate_heatmap2(mat[cl1,cl1], level + 1),
#                 branch2=generate_heatmap2(mat[cl2,cl2], level + 1)))
# }

# L <- generate_heatmap2(conmat, 0)

library(UpSetR)
library(openxlsx)
library(VennDiagram)
library(Vennerable)

get_genes <- function (tfnames, key, fract = 0.75) {

  genes_in_cluster <- lapply(tfnames, function(tf) genelist[[tf]][[key]])
  all_genes <- Reduce(union, genes_in_cluster)

  mat <- matrix(0, length(all_genes), length(tfnames), dimnames=list(all_genes, tfnames))
  for (i in seq(length(genes_in_cluster))) {
    g <- genes_in_cluster[[i]]
    mat[,i] <- as.numeric(rownames(mat) %in% g)
  }

  thr <- floor(ncol(mat) * fract)
  rs <- rowSums(mat)
  k <- rs > thr
  rownames(mat)[k]
}

find_intersections <- function(key, AR_cluster, target_clusters, outpath, width = 7) {


  AR_cluster_tfnames <- rownames(readRDS(AR_cluster))
  target_clusters_tfnames <- lapply(target_clusters, function(p) rownames(readRDS(p)))


  ar_genes <- get_genes(AR_cluster_tfnames, key)

  target_cluster_genes <- setNames(vector("list", length(target_clusters_tfnames)), names(target_clusters_tfnames))
  for(i in seq_along(target_clusters_tfnames)) {
    cluster <- target_clusters_tfnames[[i]]
    cluster_name <- names(target_clusters_tfnames)[i]

    # clusters_genes <- lapply(cluster, get_tf_genes)
    # tf_genes <- Reduce(union, clusters_genes)
    tf_genes <- get_genes(cluster, key)

    tryCatch({
      target_cluster_genes[[cluster_name]] <- tf_genes
    }, error = function(er){
        print(cluster_name)
        print(tf_genes)
        stop(er)
    })

  }

  target_cluster_genes[["AR_cluster"]] <- ar_genes

  df <- data.frame(gene=Reduce(union, target_cluster_genes))
  for(i in seq_along(target_cluster_genes)) {
    column <- names(target_cluster_genes)[i]
    df[,column] <- as.numeric(df$gene %in% target_cluster_genes[[i]])
  }

  plt <- upset(df,
    nsets = ncol(df) - 1,
    nintersects = NA,
    order.by = "freq",
    set_size.show = TRUE)

  pdf(outpath, onefile=FALSE, width = width)
  print(plt)
  dev.off()

  return(df)
}

convert_df_to_list <- function(df) {
    l <- setNames(vector("list", ncol(df) - 1), colnames(df)[2:ncol(df)])
    for (i in seq(2,ncol(df))) {
        col <- colnames(df)[i]
        l[[col]] <- df$gene[which(df[,col] == 1)]
    }
    return(l)
}

extract_intersections <- function(df, thr = 10) {
    require(Vennerable)
    l <- convert_df_to_list(df)
    v <- Venn(l)
    set_names <- names(l)
    partitions <- v@IntersectionSets
    part_lengths <- sapply(partitions, length)
    partitions <- partitions[part_lengths >= thr]
    part_names <- names(partitions)
    part_names_split <- strsplit(part_names, "")
    part_names_split <- sapply(part_names_split, function(x) {
                               x <- as.numeric(x);
                               names(x) <- set_names;
                               x <- x[x != 0]
                               n <- paste(names(x), collapse = " & ")
                               return(n) })
    index <- data.frame(sheet_name=names(partitions), clusters=part_names_split)
    partitions <- lapply(partitions, function(x) {
                             x <- as.data.frame(x)
                             x <- merge(x, mirna_map, by.x="x", by.y="mirna_name", all.x=TRUE)
                             x$map_id <- x[,1]
                             x[!is.na(x[,2]), "map_id"] <- x[!is.na(x[,2]), 2]
                             x <- merge(x, ann.gene, by.x="map_id", by.y="gene_id", all.x=TRUE)
                             x <- subset(x, select = c("x", "gene_name", "map_id"))
                             return(x)})
    partitions <- c(index=list(index), partitions)
    return(partitions)
}


# panther_enrichment <- function(geneList,
#   organism="Homo sapiens",
#   enrichmentType="process",
#   test_type="FISHER",
#   correction="FDR",
#   type="enrichment") {
#   # API doc: http://pantherdb.org/help/PANTHERhelp.jsp#V.E.

#   if(!require(httr)) install.packages("httr")
#   if(!require(data.table)) install.packages("data.table")

#   url <- "http://www.pantherdb.org/webservices/garuda/tools/enrichment/VER_2/enrichment.jsp? "

#   write(geneList, "/tmp/geneList.txt")

#   response <- POST(url,
#     body = list(organism=organism,
#       geneList=upload_file("/tmp/geneList.txt"),
#       enrichmentType=enrichmentType,
#       test_type=test_type,
#       correction=correction,
#       type=type),
#     encode = "multipart")

#   df <- fread(text=content(response, as="parsed"))
#   Sys.sleep(1)

#   return(df)

# }

processGoTable <- function(df) {

  if(!require(data.table)) install.packages(data.table)
  if (!is.data.table(df)) df <- as.data.table(df)

  dff <- df[FDR < 0.05]
  dff <- dff[,.(Genes=paste(GeneId, collapse=", ")), by=c("Id", "Name", "FDR")]
  dff <- dff[order(FDR)]
  return(dff)
}

merge_gene_list <- function(x, lfc_table, de_vect) {
  if (length(x) == 0) return(x)

  x <- as.data.frame(x)
  x <- merge(x, mirna_map, by.x="x", by.y="mirna_name", all.x=TRUE)
  x$map_id <- x[,1]
  x[!is.na(x[,2]), "map_id"] <- x[!is.na(x[,2]), 2]
  x <- merge(x, ann.gene, by.x="map_id", by.y="gene_id", all.x=TRUE)
  x <- subset(x, select = c("x", "gene_name", "map_id"))

  x$tmp <- x$x
  x$tmp <- gsub("miR", "mir", x$tmp)
  x <- merge(x, lfc_table, by.x="tmp", by.y=0, all.x=T)
  x$tmp <- NULL

  x$de <- gsub("miR", "mir", x$x) %in% de_vect

  return(x)
}


## pc_bph
root <- "/bmt-data/genomics/projects/atac_workdir/gene_expression_correlation/tf_analysis_expressed_tf_only/pc_bph"
AR_cluster <- file.path(root, "0B/conmat1B.rds")
target_clusters <- file.path(root,
  c(
    HOXB13_cluster = "0A/1B/2A/conmat3A.rds",
    CTCF_cluster   = "0A/1B/2B/conmat3B.rds",
    NKX3_cluster   = "0A/1A/2B/3B/4A/conmat5A.rds",
    TFAP4_cluster  = "0A/1A/2B/3B/4B/conmat5B.rds",
    SOX13_cluster  = "0A/1A/2B/3A/4B/conmat5B.rds",
    KMT2A_cluster  = "0A/1A/2B/3A/4A/5A/conmat6A.rds",
    TAF1_cluster   = "0A/1A/2B/3A/4A/5B/conmat6B.rds",
    MLX_cluster    = "0A/1A/2A/3B/4A/5A/conmat6A.rds"
  ))
names(target_clusters) <-   c(
    "HOXB13_cluster",
    "CTCF_cluster",
    "NKX3_cluster",
    "TFAP4_cluster",
    "SOX13_cluster",
    "KMT2A_cluster",
    "TAF1_cluster",
    "MLX_cluster")
outpath <- file.path(root,  "clusters_comparison_075.pdf")
df <- find_intersections("pc_bph", AR_cluster, target_clusters, outpath, width=21)
intersections <- extract_intersections(df)
intersections_path <- file.path(root, "genes_in_intersections.xlsx")
write.xlsx(intersections, file = intersections_path)

# compare AR_cluster and CTCF_cluster, the closest one in the tree
target_clusters <- file.path(root, c("0A/1B/conmat2B.rds", "0B/conmat1B.rds"))
names(target_clusters) <- c("other","AR_cluster")
genes_clusters_level1 <- lapply(target_clusters, function(cluster) get_genes(rownames(readRDS(cluster)), key = "pc_bph", fract = 0.3))
# 122 280

lfc_genes_pc_bph <- subset(rna_PCvBPH, select="log2FoldChange")
colnames(lfc_genes_pc_bph) <- "lfc"
lfc_mirna_pc_bph <- subset(srna_lfc, select = "log2FoldChange_pc_bph")
colnames(lfc_mirna_pc_bph) <- "lfc"
lfc_pc_bph <- rbind(lfc_genes_pc_bph,lfc_mirna_pc_bph)

de_vect_pc_bph <- c(rownames(de_PCvBPH), rownames(de_srna_PCvBPH))

genes_clusters_level1_mapped <- lapply(genes_clusters_level1, merge_gene_list, lfc_table = lfc_pc_bph, de_vect =de_vect_pc_bph )
genes_clusters_level1_mapped_path <- file.path(root, "comparison_AR_sibling.xlsx")
write.xlsx(genes_clusters_level1_mapped, file = genes_clusters_level1_mapped_path)

venn.diagram(genes_clusters_level1, filename=file.path(root, "venn_AR_sibling.png"), imagetype="png")

# enrichment_tables <- lapply(genes_clusters_level1, panther_enrichment)
# enrichment_tables <- lapply(enrichment_tables, processGoTable)


multithr_venn <- function (target_clusters, key, root) {
  thrs <- seq(0.1,1,by=0.1)
  m <- matrix(0, 4, length(thrs), dimnames=list(NULL, as.character(thrs)))
  for (thr in thrs) {

    genes_clusters_level1 <- lapply(target_clusters,
      function(cluster) get_genes(rownames(readRDS(cluster)), key = key, fract = thr))

    v <- Venn(genes_clusters_level1)
    partitions <- v@IntersectionSets

    sizes <- sapply(partitions, length)
    m[,as.character(thr)] <- sizes

  }
  rownames(m) <- names(sizes)
  m <- m[-1,]

  legend <- strsplit(rownames(m), "")
  legend <- lapply(legend, function(x) {
    x <- as.numeric(x)
    n <- names(genes_clusters_level1)[x == 1]
    paste(n, collapse = " & ")
  })

  png(file.path(root, "barplot_intersections_size.png"))
  bc <- barplot(m, col = terrain.colors(nrow(m)), main = key, xlab="Fraction of TF sharing a gene", ylab="Number of shared genes")
  H <- apply(m, 2, cumsum)
  H <- H - m / 2
  text(rep(bc, each = nrow(H)), H, labels = m)
  legend("topright", legend=legend, fill = terrain.colors(nrow(m)))
  dev.off()
}

multithr_venn(target_clusters, "pc_bph", root)



## crpc_pc
root <- "/bmt-data/genomics/projects/atac_workdir/gene_expression_correlation/tf_analysis_expressed_tf_only/crpc_pc"
AR_cluster <- file.path(root, "0B/1A/conmat2A.rds")
target_clusters <- file.path(root, c(
                                     MYC_cluster    = "0B/1A/2B/3A/4A/5A/conmat6A.rds",
                                     YY1_cluster    = "0B/1B/2B/3A/conmat4A.rds",
                                     YY1_cluster2   = "0B/1B/2B/3A/4A/conmat5A.rds",
                                     EGR3_cluster   = "0B/1B/2A/3B/conmat4B.rds",
                                     ATF3_cluster   = "0B/1B/2A/3A/conmat4A.rds",
                                     ARNT_cluster   = "0A/1B/2B/3A/4A/5A/conmat6A.rds",
                                     NCOA_cluster   = "0A/1B/2B/3A/4A/5B/conmat6B.rds",
                                     HMG20A_cluster = "0A/1B/2A/3A/4B/conmat5B.rds",
                                     AHR_cluster    = "0A/1B/2A/3A/4A/conmat5A.rds",
                                     ATF1_cluster   = "0A/1B/2A/3B/conmat4B.rds",
                                     KDM5D_cluster  = "0A/1A/2B/3A/4A/conmat5A.rds",
                                     NONO_cluster   = "0A/1A/2B/3A/4B/conmat5B.rds",
                                     ZNF281_cluster = "0A/1A/2B/3B/4A/5A/conmat6A.rds",
                                     RERE_cluster   = "0A/1A/2B/3B/4A/5B/conmat6B.rds"

))
names(target_clusters) <- c("MYC_cluster",
                            "YY1_cluster",
                            "YY1_cluster2",
                            "EGR3_cluster",
                            "ATF3_cluster",
                            "ARNT_cluster",
                            "NCOA_cluster",
                            "HMG20A_cluster",
                            "AHR_cluster",
                            "ATF1_cluster",
                            "KDM5D_cluster",
                            "NONO_cluster",
                            "ZNF281_cluster",
                            "RERE_cluster")
outpath <- file.path(root,  "clusters_comparison_075.pdf")
df <- find_intersections("crpc_pc", AR_cluster, target_clusters, outpath, width=28)
intersections <- extract_intersections(df)
intersections_path <- file.path(root, "genes_in_intersections.xlsx")
write.xlsx(intersections, file = intersections_path)

# compare AR_cluster and E2F1_cluster, the closest one in the tree
target_clusters <- file.path(root, c(sibling="0B/1B/conmat2B.rds", AR_cluster="0B/1A/conmat2A.rds"))
names(target_clusters) <- c("other","AR_cluster")
genes_clusters_level1 <- lapply(target_clusters, function(cluster) get_genes(rownames(readRDS(cluster)), key = "crpc_pc", fract = 0.3))
# 122 280

lfc_genes_crpc_pc <- subset(rna_CRPCvPC, select="log2FoldChange")
colnames(lfc_genes_crpc_pc) <- "lfc"
lfc_mirna_crpc_pc <- subset(srna_lfc, select = "log2FoldChange_crpc_pc")
colnames(lfc_mirna_crpc_pc) <- "lfc"
lfc_crpc_pc <- rbind(lfc_genes_crpc_pc,lfc_mirna_crpc_pc)

de_vect_crpc_pc <- c(rownames(de_CRPCvPC), rownames(de_srna_CRPCvPC))

genes_clusters_level1_mapped <- lapply(genes_clusters_level1, merge_gene_list, lfc_table = lfc_crpc_pc, de_vect =de_vect_crpc_pc )
genes_clusters_level1_mapped_path <- file.path(root, "comparison_AR_sibling.xlsx")
write.xlsx(genes_clusters_level1_mapped, file = genes_clusters_level1_mapped_path)

venn.diagram(genes_clusters_level1, filename=file.path(root, "venn_AR_sibling.png"), imagetype="png")
multithr_venn(target_clusters, "crpc_pc", root)

# enrichment_tables <- lapply(genes_clusters_level1, panther_enrichment)
# enrichment_tables <- lapply(enrichment_tables, processGoTable)


map_names <- function(partitions_names, n) {
  partitions_names <- strsplit(partitions_names, "")
  sapply(partitions_names, function(part) {
    part <- as.numeric(part)
    ns <- n[part == 1]
    paste(ns, collapse = " & ")
  })
}

# split the AR_cluster manually
ar_cluster_conmat <- readRDS(AR_cluster)
ar_cluster_target_tf <- c("AR", "ESR1", "FOXA1", "ERG")
# ar_cluster_target_conmat <- ar_cluster_conmat[ar_cluster_target_tf, ar_cluster_target_tf]
ar_cluster_target_genes <- get_genes(ar_cluster_target_tf, "crpc_pc", 0.3)

# take remainder of AR cluster
ar_cluster_other_tf <- rownames(ar_cluster_conmat)[!rownames(ar_cluster_conmat) %in% ar_cluster_target_tf]
# ar_cluster_other_conmat <- ar_cluster_conmat[ar_cluster_other_tf, ar_cluster_other_tf]
ar_cluster_other_genes <- get_genes(ar_cluster_other_tf, "crpc_pc", 0.3)

# take the KMT2A+TP63 cluster
kmt2a_tp63_cluster <- rownames(readRDS(file.path(root, "0B/1B/2B/conmat3B.rds")))
kmt2a_tp63_cluster_genes <- get_genes(kmt2a_tp63_cluster, "crpc_pc", 0.3)

L <- list(AR_small=ar_cluster_target_genes,
  AR_others=ar_cluster_other_genes,
  KMT2A_TP63_cluster=kmt2a_tp63_cluster_genes)
v <- Venn(L)

png(file.path(root, "venn_ArSmallcluster_vs_RemainderArCluster_vs_Kmt2aCluster.png"))
plot(v)
dev.off()

xlsx <- v@IntersectionSets[-1]
xlsx <- lapply(xlsx, merge_gene_list, lfc_table = lfc_crpc_pc, de_vect = de_vect_crpc_pc)

index <- data.frame(sheet_name=names(xlsx))
index$partition <- map_names(names(xlsx), names(L))
xlsx <- c(list(index=index), xlsx)
write.xlsx(xlsx, file=file.path(root, "ArSmallcluster_vs_RemainderArCluster_vs_Kmt2aCluster.xlsx"))




# analyze single genes
ar_foxa1_genes <- get_genes(c("AR", "FOXA1"), "crpc_pc", 0.3)
sp1_fli1_genes <- get_genes(c("SP1", "FLI1"), "crpc_pc", 0.3)
kmt2a_tp63_genes <- get_genes(c("KMT2A", "TP63"), "crpc_pc", 0.3)


L <- list(AR_FOXA1=ar_foxa1_genes,
  SP1_FLI1=sp1_fli1_genes,
  KMT2A_TP63=kmt2a_tp63_genes)
v <- Venn(L)

png(file.path(root, "venn_ArFoxa1_vs_Sp1Fli1_vs_Kmt2aTp63.png"))
plot(v)
dev.off()

xlsx <- v@IntersectionSets[-1]
xlsx <- lapply(xlsx, merge_gene_list, lfc_table = lfc_crpc_pc, de_vect = de_vect_crpc_pc)

index <- data.frame(sheet_name=names(xlsx))
index$partition <- map_names(names(xlsx), names(L))
xlsx <- c(list(index=index), xlsx)
write.xlsx(xlsx, file=file.path(root, "ArFoxa1_vs_Sp1Fli1_vs_Kmt2aTp63.xlsx"))




## peaks
root <- "/bmt-data/genomics/projects/atac_workdir/gene_expression_correlation/tf_analysis_expressed_tf_only/peaks"
AR_cluster <- file.path(root, "0B/1B/2A/conmat3A.rds")
target_clusters <- file.path(root,c(MYC_cluster     = "0B/1B/2A/3B/4A/5A/conmat6A.rds",
                                    SP2_cluster     = "0B/1A/2A/3A/4A/conmat5A.rds",
                                    BHLHE40_cluster = "0B/1A/2A/3A/4B/conmat5B.rds",
                                    HOXB13_cluster  = "0B/1A/2A/3B/conmat4B.rds",
                                    SMAD4_cluster   = "0B/1A/2B/conmat3B.rds",
                                    SMARCE1_cluster = "0A/1B/2A/3A/4A/conmat5A.rds",
                                    OSR2_cluster    = "0A/1B/2A/3A/4B/conmat5B.rds",
                                    THRB_cluster    = "0A/1B/2A/3B/4A/conmat5A.rds",
                                    SMAD5_cluster   = "0A/1B/2B/3A/4A/conmat5A.rds",
                                    HBP1_cluster    = "0A/1B/2B/3A/4B/conmat5B.rds",
                                    ZSCAN9_cluster  = "0A/1B/2B/3B/4B/conmat5B.rds",
                                    RUNX3_cluster   = "0A/1B/2B/3B/4A/5A/conmat6A.rds",
                                    HMG20A_cluster  = "0A/1B/2B/3B/4A/5B/conmat6B.rds"))
names(target_clusters) <- c("MYC_cluster",
                            "SP2_cluster",
                            "BHLHE40_cluster",
                            "HOXB13_cluster",
                            "SMAD4_cluster",
                            "SMARCE1_cluster",
                            "OSR2_cluster",
                            "THRB_cluster",
                            "SMAD5_cluster",
                            "HBP1_cluster",
                            "ZSCAN9_cluster",
                            "RUNX3_cluster",
                            "HMG20A_cluster" )
outpath <- file.path(root,  "clusters_comparison_075.pdf")
df <- find_intersections("peaks", AR_cluster, target_clusters, outpath, width = 50)
intersections <- extract_intersections(df, thr = 20)
intersections_path <- file.path(root, "genes_in_intersections.xlsx")
write.xlsx(intersections, file = intersections_path)

