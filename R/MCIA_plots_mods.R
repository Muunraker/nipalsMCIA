#' projection_plot
#'
#' @description Function to generate a projection plot MCIA results.
#' 
#' @details Plotting function for a project plot. 
#' 
#' @param mcia_result MCIA results object returned from `nipals_multiblock`
#' @param plotType Type of plot, with the following options \itemize{
#' \item `projection` - scatter plot of two orders of global and block scores (aka factors).
#' \item `projection_global` - scatter plot of two orders of global scores only (aka factors).
#' \item `gs_eigvals` - scree plot of global score eigenvalues.
#' \item `block_weights_heatmap` - a heatmap of block weightings for each order of global score.
#' }
#' @param orders Option to selecct orders of factors to plot against each other (for projection plots)
#' @param coloring Option to plot clusters/colors in projection plots, with two options:\itemize{
#' \item `none` (Default) - no clusters/colors plotted.
#' \item A character string of the column name of the `mcia_result$metadata` dataframe 
#' determining which color groupings to use (projection plots only) 
#' }
#' @param labelColors Option to specify  colors for labels.
#' @param legend_loc Option for legend location, or "none" for no legend.
#' 
#' @examples
#' data(NCI60)
#' mcia_res <- nipals_multiblock(data_blocks, metadata = metadata_NCI60, num_PCs = 10, 
#'                               plots = 'none', tol=1e-12)
#' clus_colors <- list("red", "green","blue")
#' MCIA_plots(mcia_res,'projection',orders = c(1,2), coloring="cancerType",
#'            labelColors = clus_colors, legend_loc = "bottomleft")
#' @return Displays the desired plots
#' @export
projection_plot <- function(mcia_result, plotType,
                           orders=c(1,2),
                           coloring = NULL,
                           legend_loc = "bottomleft"){
    
    ### Identifying the membership of samples within
    ### the clusters/categories of the coloring column
    # case i) no, clusters/categories were not specified 
    if(is.null(coloring)){
        clust_indexes = list(1:dim(mcia_result$global_score)[[1]])
        
    # case ii) yes, clusters/categories were specified 
    } else if(is.character(coloring)){
        
        # locating the coloring column within metadata
        col_idx <- grep(coloring, names(mcia_result$metadata))
        if(any(length(col_idx) < 1)){
            stop("Column name for coloring not found in metadata.")
        }
        
        # catching if two columns happen to have the same name.
        if (length(col_idx) > 1){
            msg = paste0('Metadata has duplicate columns for ', coloring,
                         '. Selecting the first one for plotting.')
            warning(msg)
        }
        col_idx <- col_idx[[1]] 
        
        # creating a list of samples indexes to color each cluster/category
        clust_indexes <- list()
        uniq_clusts <- unique(unlist(mcia_result$metadata[col_idx]))
        for(clust in uniq_clusts){
            clust_indexes <- c(clust_indexes, list(grep(clust, mcia_result$metadata[[col_idx]])))
        }
        names(clust_indexes) <- uniq_clusts
        
    # case iii) catch-all, coloring argument not recognized
    } else{
        message("coloring option not recognized, defaulting to black/white plotting.")
        clust_indexes = list(1:dim(mcia_result$global_score)[[1]])
    }
    
    ### Resolving the cluster colors
    plot_colors = get_metadata_colors(mcia_result, coloring = coloring)
    if (is.null(coloring)){
      plot_colors = list("black")
    }
    
    ### Plot 1 - projection plot
    if(tolower(plotType) == "projection"){
    
        # Normalize global scores to unit variance
        print("# Normalize global scores to unit variance")
        gs_norms <- apply(mcia_result$global_scores, 2, function(x){sqrt(var(x))})
        gs_normed <- t(t(mcia_result$global_scores) / gs_norms)
        gl_normed <- t(t(mcia_result$global_loadings) / gs_norms)
        gw_normed <- t(t(mcia_result$block_score_weights) / gs_norms)
        
        # Normalize block scores to unit variance
        print("# Normalize block scores to unit variance")
        bs_normed <- list()
        bl_normed <- list() 
        for(i in 1:length(mcia_result$block_scores)){
          bs_norms <-apply(mcia_result$block_scores[[i]], 2, function(x){sqrt(var(x))})
          bs_normed[[i]] <- t(t(mcia_result$block_scores[[i]]) / bs_norms)
          bl_normed[[i]] <- t(t(mcia_result$block_loadings[[i]]) / bs_norms)
        }
        
        # Getting bounds for projection plot
        print("# Getting bounds for projection plot")
        min_bs1 <- min(sapply(lapply(bs_normed, `[`,,orders[[1]]) , min)) # minimum 1st block score
        min_bs2 <- min(sapply(lapply(bs_normed, `[`,,orders[[2]]) , min)) # minimum 2nd block score
        max_bs1 <- max(sapply(lapply(bs_normed, `[`,,orders[[1]]) , max)) # maximum 1st block score
        max_bs2 <- max(sapply(lapply(bs_normed, `[`,,orders[[2]]) , max)) # maximum 2nd block score
        
        min_x <- min(c(min_bs1, min(gs_normed[,orders[[1]]]))) # minimum x coordinate in plot
        min_y <- min(c(min_bs2, min(gs_normed[,orders[[2]]]))) # minimum y coordinate in plot
        max_x <- max(c(max_bs1, max(gs_normed[,orders[[1]]]))) # maximum x coordinate in plot
        max_y <- max(c(max_bs2, max(gs_normed[,orders[[2]]]))) # maximum y coordinate in plot
        
        # Cluster 1
        print("# Cluster 1")
        sample_indexes <- clust_indexes[[1]]
        
        # Plotting global scores
        print("# Plotting global scores")
        plot(gs_normed[sample_indexes, orders[[1]]],
             gs_normed[sample_indexes, orders[[2]]],
             main = "Factor Plot",  
             xlab=paste('Factor ', orders[[1]]),
             ylab=paste('Factor ', orders[[2]]),
             col=plot_colors[[1]],
             xlim=c(min_x, max_x),
             ylim=c(min_y, max_y),
             cex = .5,pch = 16)
        grid()
        
        # Plotting block scores (shapes correspond to different blocks)
        print("# Plotting block scores (shapes correspond to different blocks)")
        for(j in 1:length(bs_normed)){
          bs_j <- bs_normed[[j]]
          points(bs_j[sample_indexes, orders[[1]]],
                 bs_j[sample_indexes, orders[[2]]], 
                 col=plot_colors[[1]],
                 cex = 1, pch = j-1)
          
          # Line segments joining block scores to central global score:
          segments(bs_j[sample_indexes, orders[[1]]],
                   bs_j[sample_indexes, orders[[2]]],
                   gs_normed[sample_indexes, orders[[1]]],
                   gs_normed[sample_indexes, orders[[2]]],
                   col=plot_colors[[1]])
        }
        
        # Cluster 2+
        if(length(clust_indexes)>1){
          for(i in 2:length(clust_indexes)){
            sample_indexes <- clust_indexes[[i]]
            points(gs_normed[sample_indexes, orders[[1]]],
                   gs_normed[sample_indexes, orders[[2]]],
                   col=plot_colors[[i]],
                   cex = .5,
                   pch = 16)
            
            for(j in 1:length(bs_normed)){
              bs_j <- bs_normed[[j]]
              points(bs_j[sample_indexes, orders[[1]]],
                     bs_j[sample_indexes, orders[[2]]],
                     col=plot_colors[[i]],
                     cex = 1,pch = j-1)
              segments(bs_j[sample_indexes, orders[[1]]],
                       bs_j[sample_indexes, orders[[2]]],
                       gs_normed[sample_indexes, orders[[1]]],
                       gs_normed[sample_indexes, orders[[2]]],
                       col=plot_colors[[i]])
            }  
          }
        }
        
        # Adding legend
        if(! tolower(legend_loc) == "none"){
            # plotting legend without clusters/categories
            if(length(plot_colors) == 1){
              print("# plotting legend without clusters/categories")
              legend(legend_loc, 
                     legend = c(names(mcia_result$block_loadings)),
                     pch = 1:length(mcia_result$block_loadings),
                     cex = 1)
            # plotting legend for clusters/categories
            } else{
              print("# plotting legend for clusters/categories")
              leg_labels = c(names(mcia_result$block_loadings),
                             names(plot_colors))
              leg_shapes = c(1:length(mcia_result$block_loadings),
                             rep(16, length(plot_colors)))
              leg_colors = c(rep("black", length(mcia_result$block_loadings)),
                             unname(unlist(plot_colors)))
              print('leg_colors:')
              print(leg_colors)
              legend(legend_loc,
                     legend = leg_labels,
                     pch = leg_shapes,
                     col = leg_colors,
                     cex = 1)
            }
        }
         
    } else if(tolower(plotType) == 'projection_global'){
    
        ### Plot 2 - projection plot global score only      
        # Normalize global scores to unit variance
        gs_norms <- apply(mcia_result$global_scores,2,function(x){sqrt(var(x))})
        gs_normed <- t(t(mcia_result$global_scores) / gs_norms)
        gl_normed <- t(t(mcia_result$global_loadings) / gs_norms)
        gw_normed <- t(t(mcia_result$block_score_weights) / gs_norms)
        
        
        # Getting bounds for projection plot
        min_x <-  min(gs_normed[,orders[[1]]]) # minimum x coordinate in plot
        min_y <-  min(gs_normed[,orders[[2]]]) # minimum y coordinate in plot
        max_x <-  max(gs_normed[,orders[[1]]]) # maximum x coordinate in plot
        max_y <-  max(gs_normed[,orders[[2]]]) # maximum y coordinate in plot
        
        sample_indexes <- clust_indexes[[1]]
        
        plot(gs_normed[sample_indexes, orders[[1]]],
             gs_normed[sample_indexes, orders[[2]]],
             main = "Global Factor Plot",  
             xlab=paste('Factor ', orders[[1]]),
             ylab=paste('Factor', orders[[2]]),
             col=plot_colors[[1]],
             xlim=c(min_x, max_x),
             ylim=c(min_y, max_y),
             pch = 16,
             cex = 0.5)
        grid()
        
        # Cluster 2+
        if(length(clust_indexes)>1){
          for(i in 2:length(clust_indexes)){
            sample_indexes <- clust_indexes[[i]]
            points(gs_normed[sample_indexes,orders[[1]]],
                   gs_normed[sample_indexes,orders[[2]]],
                   col=plot_colors[[i]],
                   cex = .5,pch = 16)
          }
        }
    }
}

#' global_scores_eigenvalues_plot
#'
#' @description Function to plot eigenvalues of scores up to num_PCs
#' 
#' @details Plotting function for eigenvalues of scores up to num_PCs
#' 
#' @param mcia_result MCIA results object returned from `nipals_multiblock`
#' @param plotType Type of plot, with the following options \itemize{
#' \item `projection` - scatter plot of two orders of global and block scores (aka factors).
#' \item `projection_global` - scatter plot of two orders of global scores only (aka factors).
#' \item `gs_eigvals` - scree plot of global score eigenvalues.
#' \item `block_weights_heatmap` - a heatmap of block weightings for each order of global score.
#' }
#' @param orders Option to selecct orders of factors to plot against each other (for projection plots)
#' @param coloring Option to plot clusters/colors in projection plots, with two options:\itemize{
#' \item `none` (Default) - no clusters/colors plotted.
#' \item A character string of the column name of the `mcia_result$metadata` dataframe 
#' determining which color groupings to use (projection plots only) 
#' }
#' @param labelColors Option to specify  colors for labels.
#' @param legend_loc Option for legend location, or "none" for no legend.
#' 
#' @examples
#' data(NCI60)
#' mcia_res <- nipals_multiblock(data_blocks, metadata = metadata_NCI60, num_PCs = 10, 
#'                               plots = 'none', tol=1e-12)
#' clus_colors <- list("red", "green","blue")
#' MCIA_plots(mcia_res,'projection',orders = c(1,2), coloring="cancerType",
#'            labelColors = clus_colors, legend_loc = "bottomleft")
#' @return Displays the desired plots
#' @export
global_scores_eigenvalues_plot <- function(mcia_result){
    
  # extract eigenvalues
  barploteigs <- unlist(mcia_result$eigvals)^2
  names(barploteigs) <- 1:length(mcia_result$eigvals)
  
  # generate barplot 
  barplot(barploteigs, xlab="Factors", cex.names = 1, 
          main = "Global Factor Score Eigenvalues ")
}

#' block_weights_heatmap
#'
#' @description Function to plot heatmap of block score weights
#' 
#' @details Plotting function for heatmap of block score weights
#' 
#' @param mcia_result MCIA results object returned from `nipals_multiblock`
#' @param plotType Type of plot, with the following options \itemize{
#' \item `projection` - scatter plot of two orders of global and block scores (aka factors).
#' \item `projection_global` - scatter plot of two orders of global scores only (aka factors).
#' \item `gs_eigvals` - scree plot of global score eigenvalues.
#' \item `block_weights_heatmap` - a heatmap of block weightings for each order of global score.
#' }
#' @param orders Option to selecct orders of factors to plot against each other (for projection plots)
#' @param coloring Option to plot clusters/colors in projection plots, with two options:\itemize{
#' \item `none` (Default) - no clusters/colors plotted.
#' \item A character string of the column name of the `mcia_result$metadata` dataframe 
#' determining which color groupings to use (projection plots only) 
#' }
#' @param labelColors Option to specify  colors for labels.
#' @param legend_loc Option for legend location, or "none" for no legend.
#' 
#' @examples
#' data(NCI60)
#' mcia_res <- nipals_multiblock(data_blocks, metadata = metadata_NCI60, num_PCs = 10, 
#'                               plots = 'none', tol=1e-12)
#' clus_colors <- list("red", "green","blue")
#' MCIA_plots(mcia_res,'projection',orders = c(1,2), coloring="cancerType",
#'            labelColors = clus_colors, legend_loc = "bottomleft")
#' @return Displays the desired plots
#' @export
block_weights_heatmap <- function(mcia_result)){
    
    # plot heatmap of block score weights
    heatmap(mcia_result$block_score_weights, Rowv = NA, Colv = NA,
            cexRow= 1.5, xlab="Factors")
}
        
    

#mcia_result = mcia_results
#coloring = 'cancerType'

#plotType = 'projection_global'
#orders = c(1,2)
#projection_plot(mcia_result, plotType,
#                orders=c(1,2), coloring = NULL,
#                legend_loc = "bottomright")

#plotType = 'projection'
#projection_plot(mcia_result, plotType,
#                orders=c(1,2), coloring = NULL,
#                legend_loc = "bottomright")

#plotType = 'projection_global'
#orders = c(1,2)
#projection_plot(mcia_result, plotType,
#                orders=c(1,2), coloring = 'cancerType',
#                legend_loc = "bottomright")

#plotType = 'projection'
#projection_plot(mcia_result, plotType,
#                orders=c(1,2), coloring = 'cancerType',
#                legend_loc = "bottomright")
    

