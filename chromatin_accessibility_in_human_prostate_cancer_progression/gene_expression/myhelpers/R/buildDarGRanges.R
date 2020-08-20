#' buildDarGRanges
#' @description a function to build the GenomicRanges representing DARs coming from all three comparison 
#' @param path_pc_bph path to bed for pc_bph comparison 
#' @param path_pc_crpc path to bed for crpc_pc comparison 
#' @param  path_crpc_bph path to bed for crpc_bph comparison 
#' @return GRanges
#' @export
buildDarGRanges <- function (path_pc_bph, path_pc_crpc, path_crpc_bph) {
  pc_bph <- importGr(path_pc_bph)
  mcols(pc_bph) <- mcols(pc_bph)[,-c(1,2)]
  pc_bph$comparison_group <- "pc_bph"
  
  crpc_pc <- importGr(path_pc_crpc)
  mcols(crpc_pc) <- mcols(crpc_pc)[,-c(1,2)]
  crpc_pc$comparison_group <- "crpc_pc"
  
  crpc_bph <- importGr(path_crpc_bph)
  mcols(crpc_bph) <- mcols(crpc_bph)[,-c(1,2)]
  crpc_bph$comparison_group <- "crpc_bph"
  
  c(pc_bph, crpc_pc, crpc_bph)
  
}
