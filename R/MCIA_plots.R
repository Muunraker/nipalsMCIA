#' MCIA_plots
#'
#' @description Function to generate descriptive plots for MCIA.
#' 
#' @details Plotting function for a variety of MCIA visualization options. 
#' 
#' @param mcia_result MCIA results object returned from `nipals_multiblock`
#' @param plotType Type of plot, with the following options \itemize{
#' \item `projection` - scatter plot of two orders of global and block scores (aka factors).
#' \item `projection_global` - scatter plot of two orders of global scores only (aka factors).
#' \item `gs_eigvals` - scree plot of global score eigenvalues.
#' \item `block_weights_heatmap` - a heatmap of block weightings for each order of global score.
#' }
#' @param orders Option to selecct orders of factors to plot against each other (for projection plots)
#' @param clusters Option for list of indices of known clusters (for projection plots).
#' @param cluserColors Option to specify cluster color strings.
#' @param legend_loc Option for legend location, or "none" for no legend.
#' 
#' @examples
#' # Plotting clusters with different colors
#' data(NCI60)
#' mcia_res <- nipals_multiblock(data_blocks,num_PCs = 10, plots = 'none', tol=1e-12)
#' CNS = 1:6; LEU = 7:12; ME = 13:21;
#' clus_list <- list(CNS, LEU, ME)
#' clus_colors <- list("red", "green","blue")
#' MCIA_plots(mcia_res,'projection',orders = c(1,2),clusters = clus_list)
#' 
#' @export
MCIA_plots <- function(mcia_result,plotType,orders=c(1,2),
                       clusters=list(1:dim(mcia_result$global_score)[[1]]),
                       clusterColors = list("red","green","blue","yellow","brown",
                                            "orange","maroon","magenta","turquoise",
                                            "darkgreen","darkblue","gold","pink"),
                       legend_loc = "bottomleft"){

  if(tolower(plotType) == 'projection'){
###  Plot 1 - projection plot
    
    # Normalize global and block scores to unit variance
    gs_norms <- apply(mcia_result$global_scores,2,function(x){sqrt(var(x))})
    gs_normed <- t(t(mcia_result$global_scores) / gs_norms)
    gl_normed <- t(t(mcia_result$global_loadings) / gs_norms)
    gw_normed <- t(t(mcia_result$block_score_weights) / gs_norms)
    
    bs_normed <- list()
    bl_normed <- list() 
    for(i in 1:length(mcia_result$block_scores)){
      bs_norms <-apply(mcia_result$block_scores[[i]],2,function(x){sqrt(var(x))})
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
    
    if(length(clusters) > length(clusterColors)){
      stop("Too many clusters - please specify a color for each cluster.")
    }
    
    # Cluster 1
    # Plotting global scores
    indx <- clusters[[1]]
    if(length(clusters)==1){
      clusterColors[[1]] = "black"
    }
    plot(gs_normed[indx,orders[[1]]],gs_normed[indx,orders[[2]]],main = "Score Projection Plot",  
         xlab=paste('Order ',orders[[1]],' Scores'), ylab=paste('Order ',orders[[2]],' Scores'),
         col=clusterColors[[1]],
         xlim=c(min_x, max_x),
         ylim=c(min_y, max_y),
         cex = .5,pch = 16)
    grid()
    # Plotting block scores (shapes correspond to different blocks)
    for(j in 1:length(bs_normed)){
      bs_j <- bs_normed[[j]]
      points(bs_j[indx,orders[[1]]],bs_j[indx,orders[[2]]], 
             col=clusterColors[[1]],cex = 1,pch = j-1)
      # Line segments joining block scores to central global score:
      segments(bs_j[indx,orders[[1]]],bs_j[indx,orders[[2]]],gs_normed[indx,orders[[1]]],
        gs_normed[indx,orders[[2]]], col=clusterColors[[1]])
    }
    
    # Cluster 2+
    if(length(clusters)>1){
      for(i in 2:length(clusters)){
        indx <- clusters[[i]]
        points(gs_normed[indx,orders[[1]]],gs_normed[indx,orders[[2]]], col=clusterColors[[i]],cex = .5,pch = 16)
        for(j in 1:length(bs_normed)){
          bs_j <- bs_normed[[j]]
          points(bs_j[indx,orders[[1]]],bs_j[indx,orders[[2]]], col=clusterColors[[i]],cex = 1,pch = j-1)
          segments(bs_j[indx,orders[[1]]],bs_j[indx,orders[[2]]],
                   gs_normed[indx,orders[[1]]],gs_normed[indx,orders[[2]]], col=clusterColors[[i]])
        }  
      }
    }

    # Adding legend
    if(!tolower(legend_loc)=="none"){
      legend(legend_loc,legend = c(names(mcia_result$block_loadings)),pch = 1:length(data_blocks)-1,
           cex = 1)
    }
  
  }else if(tolower(plotType) == 'projection_global'){
    
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
    
    if(length(clusters) > length(clusterColors)){
      stop("Too many clusters - please specify a color for each cluster.")
    }
    
    # Cluster 1
    # Plotting global scores
    indx <- clusters[[1]]
    if(length(clusters)==1){
      clusterColors[[1]] = "black"
    }
    plot(gs_normed[indx,orders[[1]]],gs_normed[indx,orders[[2]]],main = "Score Projection Plot",  
         xlab=paste('Order ',orders[[1]],' Scores'), ylab=paste('Order ',orders[[2]],' Scores'),
         col=clusterColors[[1]],
         xlim=c(min_x, max_x),
         ylim=c(min_y, max_y),
         cex = 1)
    grid()
    
    # Cluster 2+
    if(length(clusters)>1){
      for(i in 2:length(clusters)){
        indx <- clusters[[i]]
        points(gs_normed[indx,orders[[1]]],gs_normed[indx,orders[[2]]], col=clusterColors[[i]],cex = .5,pch = 16)
      }
    }
    
    
  }else if(tolower(plotType) == 'gs_eigvals'){
####  Plot 3 - Eigenvalues of scores up to num_PCs
    barploteigs <- unlist(mcia_result$eigvals)^2
    names(barploteigs) <- 1:length(mcia_result$eigvals)
    barplot(barploteigs, xlab="Global Score Order", cex.names = 1, 
            main = "Global Score Eigenvalues ")
    
  }else if(tolower(plotType) == 'block_weights_heatmap'){
#### Plot 4 - Heatmap of block score weights
    bsweights <- mcia_result$block_score_weights
    heatmap(bsweights, Rowv = NA, Colv = NA, cexRow= 1.5, xlab="Global Score Order")
    
  }else{
    stop('Unknown selection for plotType - please run help("MCIA_plots")')
  }
  
}