% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeCorrelationDarsInPromoters.R
\name{computeCorrelationDarsInPromoters}
\alias{computeCorrelationDarsInPromoters}
\title{computeCorrelationDarsInPromoters}
\usage{
computeCorrelationDarsInPromoters(proms_gid, dars_gr, atac_mat, gene_mat,
  mirna_mat, proteins_mat, samples_blacklist = NULL, debug = FALSE,
  randomize = FALSE, id_column = "dar_id", ncores = detectCores())
}
\arguments{
\item{proms_gid}{GRanges of promoters}

\item{dars_gr}{GRanges of ATAC-seq features}

\item{atac_mat}{data.table of atac levels. Need dar_id column}

\item{gene_mat}{matrix of gene expression}

\item{mirna_mat}{matrix of mirna expression}

\item{proteins_mat}{matrix of proteins expression}

\item{samples_blacklist}{character vector of blacklisted samples}

\item{debug}{if debug, work on a smaller subset of overlaps}

\item{randomize}{boolean, shuffle samples order before correlation?}

\item{id_column}{name of the id column (dar_id or peaks_id)}

\item{ncores}{number of cores to use}
}
\value{
data.frame
}
\description{
Compute correlation between genes, mirna, proteins and dars overlapping promoter ranges.
}
