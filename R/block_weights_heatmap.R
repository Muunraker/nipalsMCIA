#' block_weights_heatmap
#'
#' @description Function to plot heatmap of block score weights
#'
#' @details Plotting function for heatmap of block score weights
#'
#' @param mcia_results MCIA results object returned from `nipals_multiblock`
#'
#' @examples
#' data(NCI60)
#' mcia_results <- nipals_multiblock(data_blocks, metadata = metadata_NCI60,
#' num_PCs = 10, plots = "none", tol = 1e-12)
#' block_weights_heatmap(mcia_results)
#' @return Displays the heatmap of block weights
#' @export
block_weights_heatmap <- function(mcia_results) {
    bs_weights <- as.matrix(data.frame(mcia_results$block_score_weights))
    colnames(bs_weights) <- c(1:10)
    ComplexHeatmap::Heatmap(bs_weights, cluster_columns = FALSE, cluster_rows = FALSE, 
            name = "weight", column_title_side = "bottom", column_title = "Factor")
}
