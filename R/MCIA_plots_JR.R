#' Plotting a heatmap of global factors scores (sample v. factors)
#'
#' @description Plots a heatmap of MCIA global scores 
#' @param global_scores the global_scores matrix after running MCIA
#' @return the ggplot2 object
#' @export
global_scores_heatmap_ComplexHeatmap <- function(global_scores){
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

#' Plotting a heatmap of global_loadings versus features
#'
#' @description Plots a heatmap of MCIA global_loadings versus factors
#' @param global_loadings the global_loadings matrix after running MCIA
#' @param data_blocks the list of matrices used to calculate the MCIA factors
#' @param omic_name the name of the omic that should be plot, should be in data_blocks
#' @param select_features a vector of numbers to filter features 
#' @return the ggplot2 object
#' @export
global_loadings_heatmap_ComplexHeatmap <- function(global_loadings,
                                                   data_blocks,
                                                   omic_name, 
                                                   select_features=NULL){
    
        # get the index of the omics we need
        omics_index = which(names(data_blocks) == omic_name)
        
        # find at what index we should start in the global loadings
        if (omics_index == 1){
            global_start = 1 
        } else {
            i = 1
            global_start = 0
            while (i != omics_index){
                global_start = global_start + ncol(data_blocks[[i]])
                i = i + 1
            }
        }
        
        # global_end is global_start plus the number or rows in the 
        # data_blocks of omic_name
        global_end = global_start + ncol(data_blocks[[omics_index]]) - 1
        
        # extract data for current omic
        # transpose so rows are the latent factors and columns are the features
        omic_data = t(global_loadings[seq(global_start, global_end),])
        
        # rename rows and columns for plotting
        rownames(omic_data) = paste0("LF", seq(1, nrow(omic_data)))
        colnames(omic_data) = colnames(data_blocks[[omic_name]])
        
        # make a heatmap of the correlations
        #color_func = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
        coltitle = sprintf(sprintf("%s Features", omic_name))
        
        # filter features basd on select_features 
        if (!is.null(select_features)){
            omic_data = omic_data[, select_features]
        }
        
        # plot the heatmap
        p = ComplexHeatmap::Heatmap(omic_data, 
                name = "GL Score", 
                column_title = coltitle,
                row_title = "Latent Factors",
                row_names_gp = grid::gpar(fontsize = 7),
                show_column_names = T,
                show_row_names = T,
                row_names_side = "right"
        )    
      return(p)
}

#' Plotting a heatmap of global_loadings versus features
#'
#' @description Plots a heatmap of MCIA global_loadings versus factors
#' @param global_loadings the global_loadings matrix after running MCIA
#' @param data_blocks the list of matrices used to calculate the MCIA factors
#' @param omic_name the name of the omic that should be plot, should be in data_blocks
#' @param select_features a vector of numbers to filter features 
#' @return the ggplot2 object
#' @export
global_loadings_heatmap_ComplexHeatmap <- function(global_loadings,
                                                   data_blocks,
                                                   omic_name, 
                                                   select_features=NULL){
    
        # get the index of the omics we need
        omics_index = which(names(data_blocks) == omic_name)
        
        # find at what index we should start in the global loadings
        if (omics_index == 1){
            global_start = 1 
        } else {
            i = 1
            global_start = 0
            while (i != omics_index){
                global_start = global_start + ncol(data_blocks[[i]])
                i = i + 1
            }
        }
        
        # global_end is global_start plus the number or rows in the 
        # data_blocks of omic_name
        global_end = global_start + ncol(data_blocks[[omics_index]]) - 1
        
        # extract data for current omic
        # transpose so rows are the latent factors and columns are the features
        omic_data = t(global_loadings[seq(global_start, global_end),])
        
        # rename rows and columns for plotting
        rownames(omic_data) = paste0("LF", seq(1, nrow(omic_data)))
        colnames(omic_data) = colnames(data_blocks[[omic_name]])
        
        # make a heatmap of the correlations
        #color_func = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
        coltitle = sprintf(sprintf("%s Features", omic_name))
        
        # filter features basd on select_features 
        if (!is.null(select_features)){
            omic_data = omic_data[, select_features]
        }
        
        # plot the heatmap
        p = ComplexHeatmap::Heatmap(omic_data, 
                name = "GL Score", 
                column_title = coltitle,
                row_title = "Latent Factors",
                row_names_gp = grid::gpar(fontsize = 7),
                show_column_names = T,
                show_row_names = T,
                row_names_side = "right"
        )    
      return(p)
}