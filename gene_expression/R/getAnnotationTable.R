# load custom gff/gtf annotaion
getAnnotationTable <- function(annotation_table_path = "/bmt-data/genomics/reference/Homo_sapiens.GRCh38.90.chr.gtf",
                               chromosome_length_path = "/bmt-data/genomics/reference/GRCh38p11_chr_lengths.tsv",
                               annotation_style = "UCSC", assembly = "hg38", 
                               include_chromosomes = FALSE, include_seqlengths = TRUE) {
  suppressPackageStartupMessages(require("rtracklayer"))

  # import tables
  ann <- import.gff(annotation_table_path)
  # fix chromosome notation to chr*
  seqlevelsStyle(ann) <- annotation_style

  # set assembly
  genome(ann) <- assembly
  # ann <- sortSeqLevels(ann)

  # set chromosomes length
  if (include_seqlengths) {
    chrlength <- read.table(chromosome_length_path, header = T)
    seqlengths(ann)[match(names(seqlengths(ann)), chrlength$seqnames, nomatch = 0)] <- chrlength$length
  }

  #ann <- ann[ann$gene_name %in% rnaseq_genenames, ]
  if (!include_chromosomes){
    ann <- ann[ann$type != "chromosome", ]
  }

  return(ann)
}
