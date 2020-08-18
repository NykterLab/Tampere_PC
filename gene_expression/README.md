# Gene Expression

In the `notebooks` folder there are rmarkdown notebooks and R scripts used to: call differentially expressed genes, run the correlation analysis with ATAC-seq, detect transcription factors (TFs) binding sites from ChIP-seq data (GTRD), build and analyse the TF-gene expression network.

In the `myhelpers` folder there is are helper function written for this project. They are implemented as an R package which can be installed with the function given below. The `devtools` and `roxygen2` packages are required:

```r
recompile_package <- function (path) {
    require(devtools)
    require(roxygen2)
    original_wd <- getwd()
    setwd(path)
    document()
    setwd(dirname(path))
    install(basename(path))
    setwd(original_wd)
}

# assuming current working directory is this directory
recompile_package("./myhelpers")
```
