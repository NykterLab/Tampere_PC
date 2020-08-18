#' Import path as data.table. Rename columns to uniprot and ensgid and collapse by uniprot id. Generate 1 to many relationship table
#'
#' @param path path to table.
#' @param selected_columns columns to extract from path.
#' @return a data.table
#' @import data.table
#' @export
#' @examples
#' importMappingTable("some/path/to/mapping.txt", selected_columns = c(1,2))
#' importMappingTable("some/path/to/mapping.txt")
importMappingTable <- function (path, selected_columns = NULL) {

  require(data.table)

  if (!is.null(selected_columns)) {
    map_table <- fread(path, select = selected_columns)
  } else {
    map_table <- fread(path)
  }

  setnames(map_table, c("uniprot", "ensgid"))

  map_table[, { ensgid := paste(ensgid, collapse = "|") }, by = uniprot]

  return(map_table)
}
