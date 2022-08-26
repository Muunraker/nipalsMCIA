
#' Main NIPALS computation loop
#'
#' @description Applies the full adjusted NIPALS algorithm to generate block and 
#' global scores/loadings with the desired deflation method.  
#' 
#' @details Follows the NIPALS algorithm as described by Hanafi et. al. (2010).
#' For each order of scores/loadings, the vectors are computed via the `NIPALS_iter` 
#' function, then used to deflate the data matrix according to the desired deflation method. 
#' This process is repeated up to the desired maximum order of scores/loadings.
#' 
#' @param data_blocks a list of data frames, each in "sample" x "variable" format 
#' @param preprocMethod an option for the desired data pre-processing, either:\itemize{
#' \item `colprofile` (default) to transform the data into centered column profiles (see Meng et. al. 2014)
#' \item `none` for no pre-processing performed (NOT RECCOMENDED)
#' }
#' @param num_PCs the maximum order of scores/loadings
#' @param tol a number for the tolerance on the stopping criterion for NIPALS
#' @param maxIter a number for the maximum number of times NIPALS should iterate
#' @param deflationMethod an option for the desired deflation method, either: \itemize{
#' \item `block` deflation via block loadings (for MCIA, default)
#' \item `global` deflation via global scores (for CPCA)
#' }
#' @return a list containing the following: \itemize{
#' \item `global_scores` a matrix containing global scores as columns (NOT normalized to unit variance)
#' \item `global_loadings` a matrix containing global loadings as columns
#' \item `global_score_weights` a matrix of weights to express global scores as
#' a combination of block scores. Has dimensions "num_Blocks" by "num_PCs"
#' \item `eigvals` a matrix containing the eigenvalue for each computed global score. 
#' \item `block scores` a list of matrices, each contains the scores for one block
#' \item `block loadings` a list of matrices, each contains the loadings for one block (w/ unit length)
#' \item `block score weights` a matrix containing weights for each block score of each order used to construct the global scores.  
#' \item `preprocMethod` the preprocessing method used on the data.
#' }
#' @param plots an option to display varios plots of results: \itemize{
#' \item `all` displays plots of block scores, global scores, and eigenvalue scree plot
#' \item `global` displays only global score projections and eigenvalue scree plot
#' }
#' @examples 
#'  NIPALS_results <- nipals_multiblock(df_list, num_PCs = 2, tol = 1e-7, maxIter = 1000, deflationMethod = 'block')
#'  MCIA_result <- nipals_multiblock(df_list, num_PCs = 2)
#'  CPCA_result <- nipals_multiblock(df_list, num_PCs = 4,deflationMethod = 'global')
#' 
#' @export
nipals_multiblock <- function(data_blocks,preprocMethod='colprofile', num_PCs=10, tol=1e-12, max_iter = 1000, 
                              deflationMethod = 'block',plots="all"){
  num_blocks <- length(data_blocks)
  
  if(tolower(preprocMethod) == 'colprofile'){
    message("Performing centered column profile pre-processing...")
    data_blocks <- lapply(data_blocks,CCpreproc)
    message("Pre-processing completed.")
  }else{
    data_blocks <- lapply(data_blocks, as.matrix) # converting input data to matrix form
    message("No Pre-processing performed.")
  }
  
  # First NIPALS run
  message(paste("Computing order", 1 ,"scores"))
  nipals_result <- NIPALS_iter(data_blocks, tol)
  
  # Saving result
  global_scores <- nipals_result$global_scores # matrix containing global scores as columns
  global_loadings <- as.matrix(nipals_result$global_loadings) # matrix containing global loadings as columns
  block_score_weights <- nipals_result$block_score_weights# matrix containing block score weights as columns
  
  block_scores <- list() # list containing matrices of block scores
  block_loadings <- list() # list containing matrices of block loadings
  for(i in 1:num_blocks){
    block_scores[[i]] <- nipals_result$block_scores[,i]
    block_loadings[[i]] <- nipals_result$block_loadings[[i]]
  }
  
  # Computing block eigenvalue
  eigvals <- list(nipals_result$eigval);
  
  if(num_PCs>1){
    # generate scores/loadings up to number of PCs
    for(i in 2:num_PCs){
      message(paste("Computing order", i ,"scores"))
      # Deflate blocks
      if(tolower(deflationMethod) == 'block'){
        data_blocks <- mapply(deflate_block_bl, data_blocks, nipals_result$block_loadings)
      } else if(tolower(deflationMethod) == 'global'){
        data_blocks <- lapply(data_blocks, deflate_block_gs, gs=nipals_result$global_scores)
      }else{
        stop("Uknown option for deflation step - use 'block' or 'global'")
      }
      
      # Run another NIPALS iteration
      nipals_result <- NIPALS_iter(data_blocks, tol)
      
      
      # Save results
      global_scores <- cbind(global_scores, nipals_result$global_scores)
      global_loadings <- cbind(global_loadings, nipals_result$global_loadings)
      block_score_weights <- cbind(block_score_weights, nipals_result$block_score_weights)
      eigvals <- cbind(eigvals,nipals_result$eigval)
      
      for(j in 1:num_blocks){
        block_scores[[j]] <- cbind(block_scores[[j]], nipals_result$block_scores[,j])
        block_loadings[[j]] <- cbind(block_loadings[[j]], nipals_result$block_loadings[[j]])
      }
    }
  }
  
  # Formatting results
  names(block_scores) <- names(data_blocks)
  names(block_loadings) <- names(data_blocks)
  names(eigvals) <- paste("gs",1:num_PCs,sep = '')
  results_list <-list(global_scores, global_loadings, block_score_weights, 
                      block_scores, block_loadings, eigvals, tolower(preprocMethod))
  names(results_list) <- c('global_scores','global_loadings','block_score_weights',
                           'block_scores','block_loadings', 'eigvals','preprocMethod')
  
  # Plotting results
  if(tolower(plots) == 'all'){
    #### Plot 1 - first two scores as (x,y) coordinates
    
    
    # Normalize global and block scores to unit variance
    gs_norms <- apply(results_list$global_scores,2,function(x){sqrt(var(x))})
    gs_normed <- t(t(results_list$global_scores) / gs_norms)
    gl_normed <- t(t(results_list$global_loadings) / gs_norms)
    gw_normed <- t(t(results_list$block_score_weights) / gs_norms)
    
    bs_normed <- list()
    bl_normed <- list() 
    for(i in 1:length(results_list$block_scores)){
      bs_norms <-apply(results_list$block_scores[[i]],2,function(x){sqrt(var(x))})
      bs_normed[[i]] <- t(t(results_list$block_scores[[i]]) / bs_norms)
      bl_normed[[i]] <- t(t(results_list$block_loadings[[i]]) / bs_norms)
    }
    
    # Getting bounds for projection plot
    min_bs1 <- min(sapply(lapply(bs_normed, `[`,,1) , min)) # minimum 1st block score
    min_bs2 <- min(sapply(lapply(bs_normed, `[`,,2) , min)) # minimum 2nd block score
    max_bs1 <- max(sapply(lapply(bs_normed, `[`,,1) , max)) # maximum 1st block score
    max_bs2 <- max(sapply(lapply(bs_normed, `[`,,2) , max)) # maximum 2nd block score
    
    min_x <- min(c(min_bs1, min(gs_normed[,1]))) # minimum x coordinate in plot
    min_y <- min(c(min_bs2, min(gs_normed[,2]))) # minimum y coordinate in plot
    max_x <- max(c(max_bs1, max(gs_normed[,1]))) # maximum x coordinate in plot
    max_y <- max(c(max_bs2, max(gs_normed[,2]))) # maximum y coordinate in plot
    
    # Plotting first two global scores
    par(mfrow=c(1,2))
    plot(gs_normed[,1],gs_normed[,2],main = "First and Second Order Scores",  
         xlab="1st Order Scores", ylab="2nd Order Scores",
         col="black",
         xlim=c(min_x, max_x),
         ylim=c(min_y, max_y),
         cex = .5,pch = 16)
    grid()
    # Plotting block scores (shapes correspond to different blocks)
    for(j in 1:length(bs_normed)){
      bs_j <- bs_normed[[j]]
      points(bs_j[,1],bs_j[,2], col="black",cex = 1,pch = j-1)
      # Line segments joining block scores to central global score:
      segments(bs_j[,1],bs_j[,2],gs_normed[,1],gs_normed[,2], col="black")
    }
    legend("bottomleft",legend = c(names(data_blocks)),pch = 1:length(data_blocks)-1,
           cex = 1)
    
    
    ####  Plot 2 - Eigenvalues of scores up to num_PCs
    barploteigs <- unlist(eigvals)^2
    names(barploteigs) <- 1:num_PCs
    barplot(barploteigs, xlab="Global Score Order", cex.names = 1, 
            main = "Global Score Eigenvalues ")
    
    
    ####  Plot 3 - Block contributions to global scores
    # bsweights <- results_list$block_score_weights
    # heatmap(bsweights, Rowv = NA, Colv = NA, cexRow= 1.5, xlab="Global Score Order")
    # 
  }else if (tolower(plots) == 'global'){
    #### Plot 1 - first two scores as (x,y) coordinates
    
    
    # Normalize global scores to unit variance
    gs_norms <- apply(results_list$global_scores,2,function(x){sqrt(var(x))})
    gs_normed <- t(t(results_list$global_scores) / gs_norms)
    gl_normed <- t(t(results_list$global_loadings) / gs_norms)
    gw_normed <- t(t(results_list$block_score_weights) / gs_norms)
    
    
    # Getting bounds for projection plot
    
    min_x <-  min(gs_normed[,1]) # minimum x coordinate in plot
    min_y <-  min(gs_normed[,2]) # minimum y coordinate in plot
    max_x <-  max(gs_normed[,1]) # maximum x coordinate in plot
    max_y <-  max(gs_normed[,2]) # maximum y coordinate in plot
    
    # Plotting first two global scores
    par(mfrow=c(1,2))
    plot(gs_normed[,1],gs_normed[,2],main = "First and Second Order Global Scores",  
         xlab="1st Order Scores", ylab="2nd Order Scores",
         col="black",
         xlim=c(min_x, max_x),
         ylim=c(min_y, max_y),
         cex = 1)
    grid()
    
    
    ####  Plot 2 - Eigenvalues of scores up to num_PCs
    barploteigs <- unlist(eigvals)^2
    names(barploteigs) <- 1:num_PCs
    barplot(barploteigs, xlab="Global Score Order", cex.names = 1, 
            main = "Global Score Eigenvalues ")
    
    
    ####  Plot 3 - Block contributions to global scores
    # bsweights <- results_list$block_score_weights
    # heatmap(bsweights, Rowv = NA, Colv = NA, cexRow= 1.5, xlab="Global Score Order")
    # 
    
  }
  
return(results_list)
}
