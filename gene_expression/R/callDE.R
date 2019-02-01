# INIT

setwd("/bmt-data/genomics/projects/atac_workdir")
invisible(sapply(c(
    "purrr",
    "ggplot2",
    "grid",
    "gridExtra",
    "pheatmap",
    "BiocParallel",
    "rtracklayer",
    "GenomicAlignments",
    "DESeq2"
  ), library, character.only = TRUE, quietly  = TRUE))
BiocParallel::register(BiocParallel::MulticoreParam(parallel::detectCores()))

# SET CONSTANTS
P_THR   <- 0.05
LFC_THR <- 1
DIFF_THR <- 180

if(!file.exists("R_data/RNAdds.RData")){

  # LOAD ANNOTATION FILE

  source("R/getAnnotationTable.R")
  ann <- getAnnotationTable()
  chrM <- ann[seqnames(ann) == "chrM"]

  RNA_SEQ_FOLDER  <- "../tampere_pc/rna-seq"
  samples_blacklist <- c("BPH_337", "PC_4538", "PC_6102", "PC_6864", "CRPC_531")

  # LOAD CLINICAL DATA

  rna_sample_sheet <- read.csv("clinical_27_04_2012-1-rnaseq.csv", row.names="Sample.ID")
  atac_sample_sheet <- read.table("sample_sheet.csv", header = T, sep = ",", strip.white = T, row.names="Sample.name")

  rna_coldata <- data.frame(row.names = rownames(rna_sample_sheet),
                            Type = trimws(rna_sample_sheet$Tissue.Type),
                            RNA.Isolation.method = trimws(rna_sample_sheet$RNA.Isolation.protocol))
  rna_coldata$Type <- gsub(pattern = "BPH.*",
                           replacement = "BPH",
                           x = rna_coldata$Type)
  rna_coldata$RNA.Isolation.method <- gsub(pattern = "Trizol.*",
                                           replacement = "Trizol",
                                           x = rna_coldata$RNA.Isolation.method)
  rna_coldata$RNA.Isolation.method <- gsub(pattern = "\\s",
                                           replacement = "_",
                                           x = rna_coldata$RNA.Isolation.method)
  rna_coldata$Type <- factor(rna_coldata$Type)
  rna_coldata$RNA.Isolation.method <- factor(rna_coldata$RNA.Isolation.method)
  rna_coldata$sample <- rownames(rna_coldata)
  rna_coldata <- DataFrame(rna_coldata)



  # LOAD ACTUAL DATA

  rna_samples <- rownames(rna_sample_sheet)
  rna <- rna_samples %>%
          purrr::map_chr(~ file.path(getwd(), RNA_SEQ_FOLDER, paste(., "bamReadsPerGene", "out", "tab", sep = "."))) %>%
          purrr::map(safely(~ read.table(., row.names = 1, skip = 4))) %>%
          purrr::set_names(rna_samples) %>%
          purrr::transpose()
  is_ok <- rna$error %>% purrr::map_lgl(is_null)
  x <- rna$result[is_ok] %>%
        purrr::map(~ {
            v <- .[,1]
            names(v) <- rownames(.)
            v
          }) %>%
        as.data.frame()
  x <- as.matrix(x)

  rowranges <- ann[ann$gene_id %in% rownames(x) & ann$type == "gene"]

  x <- x[match(rownames(x), rowranges$gene_id), match(colnames(x), rownames(rna_coldata))]
  rna_se <- SummarizedExperiment(assays = SimpleList(counts = x), rowRanges = rowranges, colData = rna_coldata)



  # RUN DESEQ2 ANALYSIS AND GENERATE NORMALIZED DATA MATRIX

  rna_dds <- DESeqDataSet(rna_se, design = ~ RNA.Isolation.method + Type)

  rs <- rowSums(assay(rna_dds))
  # histogram!!
  lower_quartile_thr <- quantile(rs, 0.25)
  rna_dds <- rna_dds[rs > lower_quartile_thr, ] # remove lower quartile

  rna_dds <- rna_dds[,!colnames(rna_dds) %in% samples_blacklist] # remove low quality samples
  rna_dds <- rna_dds[,colnames(rna_dds) %in% rownames(atac_sample_sheet)] # keep sample for which we have atac-seq data
  rna_dds <- collapseReplicates(rna_dds, rna_dds$sample) # collapse technical replicates

  rna_dds <- DESeq(rna_dds, fitType = "local",
  #                    minReplicatesForReplace = 3,
                   parallel = F,
                   quiet = F)

  beta <- coef(rna_dds)[, "RNA.Isolation.method_Trizol_vs_Qiagen_AllPrep"]
  beta[is.na(beta)] <- 0

  p <- colData(rna_dds)$RNA.Isolation.method
  p <- sapply(p, function(each) ifelse(each == "Trizol", TRUE, FALSE))
  p <- matrix(rep(p, times = length(beta)), nrow = length(beta), byrow = T)

  sf <- sizeFactors(rna_dds)
  sf <- matrix(rep(sf, times = nrow(rna_dds)), nrow = nrow(rna_dds), byrow = T)

  norm.counts <- counts(rna_dds, normalized = F) / sf
  norm.counts <- log2(norm.counts + 1) - (beta * p)

  assay(rna_dds, "norm.counts") <- norm.counts



  # FILTER FOR PROTEIN CODING GENES LOCATED IN CHROMOSOMES OTHER THAN MITHOCONDRIAL

  rna_dds <- rna_dds[rowRanges(rna_dds)$gene_biotype == "protein_coding" &
                     !rowRanges(rna_dds)$allZero &
                     !(rowRanges(rna_dds)$gene_id %in% chrM$gene_id), ]


  save(rna_dds, file = "R_data/RNAdds.RData")
  save(x, file = "R_data/RNA_rcMatrix.RData")
} else {
  load("R_data/RNAdds.RData")
}

# GET DEG

get_deg <- function (contrast) {

  sample_set1 <- grepl(x = colnames(rna_dds), pattern = sprintf("^%s", contrast[2]))
  sample_set2 <- grepl(x = colnames(rna_dds), pattern = sprintf("^%s", contrast[3]))

  res <- results(rna_dds, contrast = contrast, alpha = P_THR)
  res$diff <- apply(2 ^ assay(rna_dds, "norm.counts")[, sample_set1], 1, median, na.rm = T) - apply(2 ^ assay(rna_dds, "norm.counts")[, sample_set2], 1, median, na.rm = T)

  sbst <- subset(res, padj < P_THR & abs(log2FoldChange) > LFC_THR & abs(diff) > DIFF_THR)

  return(list(res, sbst))
}

l <- get_deg(c("Type", "PC",   "BPH"))
rna_PCvBPH <- l[[1]]
de_PCvBPH <- l[[2]]

l <- get_deg(c("Type", "CRPC",   "PC"))
rna_CRPCvPC <- l[[1]]
de_CRPCvPC <- l[[2]]

l <- get_deg(c("Type", "CRPC",   "BPH"))
rna_CRPCvBPH <- l[[1]]
de_CRPCvBPH <- l[[2]]


save(rna_PCvBPH, de_PCvBPH,
  rna_CRPCvPC, de_CRPCvPC,
  rna_CRPCvBPH, de_CRPCvBPH, file = "R_data/RNA_differentiallyExpressed.RData")




