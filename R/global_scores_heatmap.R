#' Plotting a heatmap of global factors scores (sample v. factors)
#'
#' @description Plots a heatmap of MCIA global scores 
#' @param global_scores the global_scores matrix after running MCIA
#' @return the ggplot2 object
#' @export
global_scores_heatmap <- function(mcia_result, coloring=1){
    
    # extract color type and create a color pallette
    color_types <- mcia_result$metadata[,coloring]
    cat_colors <- scales::viridis_pal(option = "C")(length(unique(color_types)))
    names(cat_colors) <- unique(color_types)
    
    # extract the column names (if not provided via coloring
    if (typeof(coloring) == "double"){
        color_col = names(mcia_result$metadata)[coloring]
    }
    
    # add the colors to a HeatmapAnnotation
    right_colors = list(cat_colors)
    names(right_colors) = c(color_col)
    row_ha = ComplexHeatmap::HeatmapAnnotation(which = "row", 
                               ColorType = color_types,
                               col = right_colors,
                               annotation_label = color_col,
                               show_annotation_name = FALSE)
    
    # extract and re-label global scores
    global_scores = mcia_result$global_scores
    colnames(global_scores) = paste0('F', seq(1, ncol(global_scores)))
    
    p = ComplexHeatmap::Heatmap(global_scores, 
                                name = "GS Score", 
                                column_title = "Factors",
                                row_title = "Samples",
                                row_names_gp = grid::gpar(fontsize = 7),
                                show_column_names = T,
                                show_row_names = T,
                                row_names_side = "right",
                                right_annotation = row_ha
    )    
    return(p)
}