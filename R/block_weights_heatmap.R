#' block_weights_heatmap
#'
#' @description Function to plot heatmap of block score weights
#'
#' @details Plotting function for heatmap of block score weights
#'
#' @param mcia_results MCIA results object returned from `nipals_multiblock`
#' @param separate_rects a boolean indicating whether or not to add white lines between the rectangles for visual clarity
#'
#' @examples
#' data(NCI60)
#' data_blocks_mae <- simple_mae(data_blocks,row_format="sample",
#'                               colData=metadata_NCI60)
#' mcia_results <- nipals_multiblock(data_blocks_mae, num_PCs = 10,
#'                                   plots = "none", tol = 1e-12)
#' block_weights_heatmap(mcia_results)
#' @return heatmap object containing the block weights as a heatmap
#' @export
block_weights_heatmap <- function(mcia_results, separate_rects = TRUE) {
    bs_weights <- as.matrix(data.frame(mcia_results@block_score_weights))
    colnames(bs_weights) <- seq_len(ncol(bs_weights))

    # create a graphic parameters object for drawing the rectangles
    if (separate_rects) {
        rect_gp <- grid::gpar(col = "white", lwd = 1)
    } else {
        rect_gp <- grid::gpar(col = NA) # default
    }

    ComplexHeatmap::Heatmap(matrix = bs_weights,
                            name = "Weight",
                            rect_gp = rect_gp,
                            column_title = "Factors",
                            column_title_side = "bottom",
                            cluster_rows = FALSE,
                            cluster_columns = FALSE,
                            column_names_rot = 0,
                            column_names_centered = TRUE)
}
