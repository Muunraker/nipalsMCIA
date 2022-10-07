#' Assigning colors to different omics (default: color-blindness friendly)
#'
#' @description Creates a list of omics and associated colors for plotting
#' 
#' @param mcia_results object returned from nipals_multiblock() function
#' @return List of omics with assigned colors
#' @examples
#' colors_omics <- get_colors(mcia_results)
#' @importFrom scales viridis_pal
#' @export
get_colors <- function(mcia_result, color_func=scales::viridis_pal) {
  omic_list <- names(mcia_result$block_loadings)
  colors_omics <- color_func(option = "D")(length(omic_list))
  names(colors_omics) <- omic_list
  return(colors_omics)
}

#' Assigning colors to different omics (default: color-blindness friendly)
#'
#' @description Creates a list of omics and associated colors for plotting
#' 
#' @param mcia_results object returned from nipals_multiblock() function
#' @return List of omics with assigned colors
#' @examples
#' colors_omics <- get_colors(mcia_results)
#' @importFrom scales viridis_pal
#' @export
get_metadata_colors <- function(mcia_result, coloring,
                                color_func=scales::viridis_pal) {
  meta_list <- unique(mcia_result$metadata[,coloring])
  colors_meta <- color_func(option = "D")(length(meta_list))
  names(colors_meta) <- meta_list
  return(colors_meta)
}
