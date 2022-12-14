#' Plotting a heatmap of global factors scores (sample v. factors)
#'
#' @description Plots a heatmap of MCIA global scores
#' @param mcia_results the mcia object matrix after running MCIA, must also
#'   contain metadata with columns corresponding to color_col
#' @param color_col an integer or string specifying the column that will be
#'   used for color_col
#' @param color_pal a list of colors or function which returns a list of colors
#' @param color_pal_params a list of parameters for the color function
#' @return ComplexHeatmap object
#' @export
global_scores_heatmap <- function(mcia_results, color_col = NULL,
                                  color_pal = scales::viridis_pal,
                                  color_pal_params = list(option = "D")) {
    # extract the global scores
    global_scores <- mcia_results$global_scores

    # extract and re-label global scores
    colnames(global_scores) <- paste0("F", seq(1, ncol(global_scores)))

    if (!is.null(color_col)) {
        # check if metadata is within the MCIA results
        if (!"metadata" %in% names(mcia_results)) {
            msg <- paste("Cannot produce a sample color_col without metadata.",
                         "Please re-run global_scores_heatmap after setting",
                         "this information.")
            stop(msg)
        }

        # extract the column names (if not provided via color_col)
        if (typeof(color_col) == "double") {
            color_col <- names(mcia_results$metadata)[color_col]
        }

        # extract color type and create a color palette list
        cat_values <- mcia_results$metadata[, color_col]

        cat_colors <- get_metadata_colors(mcia_results,
                                          color_col = color_col,
                                          color_pal = color_pal,
                                          color_pal_params = color_pal_params)
        cat_colors_list <- list(ColorType = cat_colors)

        # add the colors list to a HeatmapAnnotation obj
        row_ha <- ComplexHeatmap::HeatmapAnnotation(which = "row",
                                                    ColorType = cat_values,
                                                    col = cat_colors_list,
                                                    annotation_label = color_col,
                                                    show_annotation_name = FALSE
            )

        p <- ComplexHeatmap::Heatmap(global_scores, name = "GS Score",
                                     column_title = "Factors",
                                     row_title = "Samples",
                                     row_names_gp = grid::gpar(fontsize = 7),
                                     show_column_names = TRUE,
                                     show_row_names = TRUE,
                                     row_names_side = "right",
                                     right_annotation = row_ha
        )
    }
    else {
        p <- ComplexHeatmap::Heatmap(global_scores, name = "GS Score",
                                     column_title = "Factors",
                                     row_title = "Samples",
                                     row_names_gp = grid::gpar(fontsize = 7),
                                     show_column_names = TRUE,
                                     show_row_names = TRUE,
                                     row_names_side = "right"
        )
    }

    return(p)
}
