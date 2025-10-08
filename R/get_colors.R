#' Assigning colors to different omics
#'
#' @description Creates a list of omics and associated colors for plotting. The
#' default palette was chosen to be color-blindness friendly.
#'
#' @param mcia_results object returned from nipals_multiblock() function
#' @param color_pal a function which returns color palettes (e.g. scales)
#' @param color_pal_params list of parameters for the corresponding function
#' @return List of omics with assigned colors
#' @examples
#' data(NCI60)
#' data_blocks_mae <- simple_mae(data_blocks, row_format = "sample",
#'                               colData = metadata_NCI60)
#' mcia_results <- nipals_multiblock(data_blocks_mae, num_PCs = 10,
#'                                  plots = "none", tol = 1e-12)
#' colors_omics <- get_colors(mcia_results)
#' @importFrom scales viridis_pal
#' @importFrom methods is
#' @export
get_colors <- function(mcia_results, color_pal = scales::viridis_pal,
                       color_pal_params = list()) {
    omic_list <- names(mcia_results@block_loadings)

    if (is(color_pal, "function")) {
        colors_omics <- do.call(color_pal, color_pal_params)(length(omic_list))

    } else if (is(color_pal, "character")) {
        colors_omics <- color_pal
    }

    names(colors_omics) <- omic_list
    return(colors_omics)
}
