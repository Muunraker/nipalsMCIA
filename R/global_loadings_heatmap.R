
#' Plotting a heatmap of global_loadings versus features
#'
#' @description Plots a heatmap of MCIA global_loadings versus factors
#' @param global_loadings the global_loadings matrix after running MCIA
#' @param data_blocks the list of matrices used to calculate the MCIA factors
#' @param omic_name the name of the omic that should be plot, should be
#' in data_blocks
#' @param select_features a vector of numbers to filter features
#' @return the ggplot2 object
#' @export
global_loadings_heatmap <- function(global_loadings, data_blocks,
                                     omic_name, select_features = NULL) {
  # get the index of the omics we need
  omics_index <- which(names(data_blocks) == omic_name)

  # find at what index we should start in the global loadings
  if (omics_index == 1) {
    global_start <- 1

  } else {
    i <- 1
    global_start <- 0
    while (i != omics_index) {
      global_start <- global_start + ncol(data_blocks[[i]])
      i <- i + 1
    }
  }

  # global_end is global_start plus the number or rows in the
  # data_blocks of omic_name
  global_end <- global_start + ncol(data_blocks[[omics_index]]) - 1

  # extract data for current omic
  # transpose so rows are the latent factors and columns are the features
  omic_data <- t(global_loadings[seq(global_start, global_end), ])

  # rename rows and columns for plotting
  rownames(omic_data) <- paste0("LF", seq(1, nrow(omic_data)))
  colnames(omic_data) <- colnames(data_blocks[[omic_name]])

  # make a heatmap of the correlations
  coltitle <- sprintf(sprintf("%s Features", omic_name))

  # filter features basd on select_features
  if (!is.null(select_features)) {
    omic_data <- omic_data[, select_features]
  }

  # plot the heatmap
  p <- ComplexHeatmap::Heatmap(matrix = omic_data,
                               name = "GL Score",
                               row_title = "Latent Factors",
                               column_title = coltitle,
                               row_names_side = "right",
                               show_row_names = TRUE,
                               row_names_gp = grid::gpar(fontsize = 7),
                               show_column_names = TRUE)

  return(p)
}
