#' Plotting a heatmap of global factors scores (sample v. factors)
#'
#' @description Plots a heatmap of MCIA global scores 
#' @param mcia_result the mcia object matrix after running MCIA, must also
#'   contain metadata with columns corresponding to coloring
#' @param coloring an integer or string specifying the column that will be
#'   used for coloring
#' @return ComplexHeatmap object
#' @export
global_scores_heatmap <- function(mcia_result, coloring=NULL){
    
    
    # extract the global scores
    global_scores = mcia_result$global_scores
    
    # extract and re-label global scores
    colnames(global_scores) = paste0('F', seq(1, ncol(global_scores)))
    
    if (! is.null(coloring)){
        
        # check if metadata is within the MCIA results
        if (! "metadata" %in% names(mcia_results)){
            msg = paste("Cannot produce a sample coloring without metadata.",
                        "Please re-run global_scores_heatmap after setting",
                        "this information.")
            stop(msg)
        }
    
        # extract the column names (if not provided via coloring)
        if (typeof(coloring) == "double"){
            coloring = names(mcia_result$metadata)[coloring]
        }
        
        # extract color type and create a color palette list
        cat_values <- mcia_result$metadata[,coloring]
        cat_colors <- get_metadata_colors(mcia_result, coloring = coloring)
        cat_colors_list <- list(ColorType = cat_colors)
        
        # add the colors list to a HeatmapAnnotation obj
        row_ha = ComplexHeatmap::HeatmapAnnotation(which = "row", 
                                   ColorType = cat_values,
                                   col = cat_colors_list,
                                   annotation_label = coloring,
                                   show_annotation_name = FALSE)
        
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
    }else{
        p = ComplexHeatmap::Heatmap(global_scores, 
                                    name = "GS Score", 
                                    column_title = "Factors",
                                    row_title = "Samples",
                                    row_names_gp = grid::gpar(fontsize = 7),
                                    show_column_names = T,
                                    show_row_names = T,
                                    row_names_side = "right"
        )    
    }
    return(p)
}