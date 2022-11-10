#' Assigning colors to different omics (default: color-blindness friendly)
#'
#' @description Creates a list of omics and associated colors for plotting
#' 
#' @param mcia_results object returned from nipals_multiblock() function
#' @param color_pal a function which returns color palettes (e.g. scales)
#' @param color_pal_params list of parameters for the corresponding function
#' @return List of omics with assigned colors
#' @examples
#' colors_omics <- get_colors(mcia_results)
#' @importFrom scales viridis_pal
#' @export
get_colors <- function(mcia_results, color_pal=scales::viridis_pal,
                                color_pal_params=list(option="D")) {
  omic_list <- names(mcia_results$block_loadings)
  
  if (is(color_pal, "function")){
      colors_omics <- do.call(color_pal, color_pal_params)(length(omic_list))
  } else if (is(color_pal, "character")) {
      colors_omics = color_pal
  }
  
  names(colors_omics) <- omic_list
  return(colors_omics)
}

#' Assigning colors to different values of a metadata column 
#' (default: color-blindness friendly)
#'
#' @description Creates a list of metadata columns and associated colors
#'   for plotting
#' @param mcia_results object returned from nipals_multiblock() function
#' @param color_col an integer or string specifying the column that will be
#'   used for color_col
#' @param color_pal a function which returns color palettes (e.g. scales)
#' @param color_pal_params list of parameters for the corresponding function
#' @return List of metadata columns with assigned colors
#' @examples
#' colors_omics <- get_metadata_colors(mcia_results, "cancerType")
#' @importFrom scales viridis_pal
#' @export
get_metadata_colors <- function(mcia_results, color_col,
                                color_pal=scales::viridis_pal,
                                color_pal_params=list(option="E")) {
    
  meta_list <- unique(mcia_results$metadata[,color_col])
  
  if (is(color_pal, "function")){
      colors_meta <- do.call(color_pal, color_pal_params)(length(meta_list))
  } else if (is(color_pal, "character")){
      colors_meta = color_pal
  }
  
  names(colors_meta) <- meta_list
  return(colors_meta)
  
}
