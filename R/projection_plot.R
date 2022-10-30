#' projection_plot
#'
#' @description Function to generate a projection plot of MCIA results.
#' 
#' @details Plotting function for a projection plot. 
#' 
#' @param mcia_results MCIA results object returned from `nipals_multiblock`
#' @param plotType Type of plot, with the following options \itemize{
#' \item `projection` - scatter plot of two orders of global and block scores (aka factors).
#' \item `projection_global` - scatter plot of two orders of global scores only (aka factors).
#' }
#' @param orders Option to select orders of factors to plot against each other (for projection plots)
#' @param color_col Option to plot clusters/colors in projection plots, with two options:\itemize{
#' \item `none` (Default) - no clusters/colors plotted.
#' \item A character string of the column name of the `mcia_result$metadata` dataframe 
#' determining which color groupings to use (projection plots only) 
#' }
#' @param color_col an integer or string specifying the column that will be
#'   used for color_col
#' @param color_pal a list of colors or function which returns a list of colors
#' @param color_pal_params a list of parameters for the color function
#' @param legend_loc Option for legend location, or "none" for no legend.
#' 
#' @examples
#' data(NCI60)
#' mcia_results <- nipals_multiblock(data_blocks, metadata = metadata_NCI60,
#'                               num_PCs = 10, plots = "none", tol=1e-12)
#' projection_plot(mcia_results, "projection", orders = c(1,2),
#'   color_col = "cancerType", legend_loc = "bottomright")
#' @return Displays the desired plots
#' @export
projection_plot <- function(mcia_result, plotType, orders=c(1,2),
                               color_col=NULL,
                               color_pal=scales::viridis_pal, 
                               color_pal_params=list(option="D"),
                               legend_loc = "bottomleft"){
    
    ### Identifying the membership of samples within
    ### the clusters/categories of the color_col column
    # case i) no color_col, clusters/categories were not specified 
    if(is.null(color_col)){
        clust_indexes = list(1:dim(mcia_result$global_score)[[1]])
        
    # case ii) yes color_col, clusters/categories were specified 
    } else if(is.character(color_col)){
        
        # locating the color_col column within metadata
        col_idx <- grep(color_col, names(mcia_result$metadata))
        if(any(length(col_idx) < 1)){
            stop("Column name for color_col not found in metadata.")
        }
        
        # catching if two columns happen to have the same name.
        if (length(col_idx) > 1){
            msg = paste0("Metadata has duplicate columns for ", color_col,
                         ". Selecting the first one for plotting.")
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
        
    # case iii) catch-all, color_col argument not recognized
    } else{
        message("color_col option not recognized, defaulting to black/white plotting.")
        clust_indexes = list(1:dim(mcia_result$global_score)[[1]])
    }
    
    ### Resolving the cluster colors
    plot_colors = get_metadata_colors(mcia_result, color_col = color_col,
                                      color_pal = color_pal,
                                      color_pal_params = color_pal_params)
    if (is.null(color_col)){
      plot_colors = list("black")
    }
    
    ### Plot 1 - projection plot
    if(tolower(plotType) == "projection"){
    
        # Normalize global scores to unit variance
        gs_norms <- apply(mcia_result$global_scores, 2, function(x){sqrt(var(x))})
        gs_normed <- t(t(mcia_result$global_scores) / gs_norms)
        gl_normed <- t(t(mcia_result$global_loadings) / gs_norms)
        gw_normed <- t(t(mcia_result$block_score_weights) / gs_norms)
        
        # Normalize block scores to unit variance
        bs_normed <- list()
        bl_normed <- list() 
        for(i in 1:length(mcia_result$block_scores)){
          bs_norms <-apply(mcia_result$block_scores[[i]], 2, function(x){sqrt(var(x))})
          bs_normed[[i]] <- t(t(mcia_result$block_scores[[i]]) / bs_norms)
          bl_normed[[i]] <- t(t(mcia_result$block_loadings[[i]]) / bs_norms)
        }
        
        # Getting bounds for projection plot
        min_bs1 <- min(sapply(lapply(bs_normed, `[`,,orders[[1]]) , min)) # minimum 1st block score
        min_bs2 <- min(sapply(lapply(bs_normed, `[`,,orders[[2]]) , min)) # minimum 2nd block score
        max_bs1 <- max(sapply(lapply(bs_normed, `[`,,orders[[1]]) , max)) # maximum 1st block score
        max_bs2 <- max(sapply(lapply(bs_normed, `[`,,orders[[2]]) , max)) # maximum 2nd block score
        
        min_x <- min(c(min_bs1, min(gs_normed[,orders[[1]]]))) # minimum x coordinate in plot
        min_y <- min(c(min_bs2, min(gs_normed[,orders[[2]]]))) # minimum y coordinate in plot
        max_x <- max(c(max_bs1, max(gs_normed[,orders[[1]]]))) # maximum x coordinate in plot
        max_y <- max(c(max_bs2, max(gs_normed[,orders[[2]]]))) # maximum y coordinate in plot
        
        # Cluster 1
        sample_indexes <- clust_indexes[[1]]
        
        # Plotting global scores
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
              legend(legend_loc, 
                     legend = c(names(mcia_result$block_loadings)),
                     pch = 0:length(mcia_result$block_loadings),
                     cex = 1)
            # plotting legend for clusters/categories
            } else{
              leg_labels = c(names(mcia_result$block_loadings),
                             names(plot_colors))
              leg_shapes = c(1:length(mcia_result$block_loadings),
                             rep(16, length(plot_colors)))-1
              leg_colors = c(rep("black", length(mcia_result$block_loadings)),
                             unname(unlist(plot_colors)))
              legend(legend_loc,
                     legend = leg_labels,
                     pch = leg_shapes,
                     col = leg_colors,
                     cex = 1)
            }
        }
         
    } else if(tolower(plotType) == "projection_global"){
    
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
             xlab=paste("Factor ", orders[[1]]),
             ylab=paste("Factor ", orders[[2]]),
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
        
        # Adding legend
        # plotting legend for clusters/categories
        if(!is.null(color_col)){
            leg_labels = c(names(plot_colors))
            leg_shapes = c(rep(16, length(plot_colors)))
            leg_colors = c(unname(unlist(plot_colors)))
            legend(legend_loc,
                   legend = leg_labels,
                   pch = leg_shapes,
                   col = leg_colors,
                   cex = 1)
        }
    }
}