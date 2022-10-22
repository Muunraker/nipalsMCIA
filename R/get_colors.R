#' Assigning colors to different omics (default: color-blindness friendly)
#'
#' @description Creates a list of omics and associated colors for plotting
#' 
#' @param mcia_result object returned from nipals_multiblock() function
#' @param color_func a function which returns color palettes (e.g. scales)
#' @param color_params list of parameters for the corresponding function
#' @return List of omics with assigned colors
#' @examples
#' colors_omics <- get_colors(mcia_results)
#' @importFrom scales viridis_pal
#' @export
get_colors <- function(mcia_result, color_func=scales::viridis_pal,
                                color_params=list(option="D")) {
  omic_list <- names(mcia_result$block_loadings)
  colors_omics <- do.call(color_func, color_params)(length(omic_list))
  names(colors_omics) <- omic_list
  return(colors_omics)
}

#' Assigning colors to different values of a metadata column 
#' (default: color-blindness friendly)
#'
#' @description Creates a list of metadata columns and associated colors
#'   for plotting
#' @param mcia_result object returned from nipals_multiblock() function
#' @param color_col an integer or string specifying the column that will be
#'   used for color_col
#' @param color_func a function which returns color palettes (e.g. scales)
#' @param color_params list of parameters for the corresponding function
#' @return List of metadata columns with assigned colors
#' @examples
#' colors_omics <- get_metadata_colors(mcia_results, "cancerType")
#' @importFrom scales viridis_pal
#' @export
get_metadata_colors <- function(mcia_result, color_col,
                                color_func=scales::viridis_pal,
                                color_params=list(option="D")) {
  meta_list <- unique(mcia_result$metadata[,color_col])
  colors_meta <- do.call(color_func, color_params)(length(meta_list))
  names(colors_meta) <- meta_list
  return(colors_meta)
}
