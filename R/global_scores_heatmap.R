#' Plotting a heatmap of global factors scores (sample v. factors)
#'
#' @description Plots a heatmap of MCIA global scores 
#' @param global_scores the global_scores matrix after running MCIA
#' @return the ggplot2 object
#' @export
global_scores_heatmap <- function(global_scores){
    colnames(global_scores) = paste0('F', seq(1, ncol(global_scores)))
    p = ComplexHeatmap::Heatmap(global_scores, 
                                name = "GS Score", 
                                column_title = "Factors",
                                row_title = "Samples",
                                row_names_gp = grid::gpar(fontsize = 7),
                                show_column_names = T,
                                show_row_names = T,
                                row_names_side = "right"
    )    
    return(p)
}