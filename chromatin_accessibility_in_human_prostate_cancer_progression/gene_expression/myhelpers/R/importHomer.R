#' importHomer
#' @description a function to load Homer annotation table as data.table
#' @param path path to the annotation table
#' @param comparison comparison name 
#' @param id_name name to assign to the first column
#' @export
importHomer <- function (path, comparison = NULL, id_name = "dar_id") {
  homer <- fread(path)
  setnames(homer, names(homer)[1], id_name, skip_absent = TRUE)
  set(homer, i = NULL, j = id_name, list( paste(homer$Chr, trimws(format(homer$Start - 1, scientific = FALSE)), homer$End, sep = "_") ))
  if (!is.null(comparison)){
    homer[,comparison_group := comparison]
  }
  return(homer)
}
