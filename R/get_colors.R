#' Assigning colors to different omics (color-blindness friendly)
#'
#' @description Creates a list of omics and associated colors for plotting
#' 
#' @param mcia_out object returned from nipals_multiblock() function
#' @return List of omics with assigned colors
#' @examples
#' colors_omics <- get_colors(mcia_out)
#' @importFrom scales viridis_pal
#' @export

get_colors <- function(mcia_out) {
  omic_list <- names(mcia_out$block_loadings)
  colors_omics <- scales::viridis_pal(option = "D")(length(omic_list))
  names(colors_omics) <- omic_list
  return(colors_omics)
}
