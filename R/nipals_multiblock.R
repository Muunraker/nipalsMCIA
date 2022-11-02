
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
#' @param max_iter a number for the maximum number of times NIPALS should iterate
#' @param metadata a data frame containing metadata (i.e. sample labels) for each sample in the dataframe.
#' May have multiple columns, but rows and row names must match the data frames in `data_blocks`.
#' @param color_col Optional argument with the column name of the `metadata` data frame used to define plotting colors
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
#' \item `none` does not display plots
#' }
#' \item `block_variances` a list of variances of each block AFTER NORMALIZATION OPTION APPLIED
#' \item `metadata` the metadata dataframe supplied wuith the `metadata` argument.
#' @examples 
#'  NIPALS_results <- nipals_multiblock(df_list, num_PCs = 10, tol = 1e-12, maxIter = 1000, 
#'                                    preprocMethod='colprofile', deflationMethod = 'block')
#'  MCIA_result <- nipals_multiblock(df_list, num_PCs = 2)
#'  CPCA_result <- nipals_multiblock(df_list, num_PCs = 4,deflationMethod = 'global')
#' 
#' @export
nipals_multiblock <- function(data_blocks,preprocMethod='colprofile', num_PCs=10, tol=1e-12, max_iter = 1000,
                              metadata = NULL, color_col = NULL, deflationMethod = 'block',plots="all"){
  num_blocks <- length(data_blocks)
  omics_names <- names(data_blocks)
  
  # Check for omics names and assign generic name if null
  if(is.null(omics_names)){
    omics_names <- paste("omic",1:length(data_blocks),sep="")
    names(data_blocks) <- omics_names
  }
  
  # Formatting feature labels to include omic type
  for(i in 1:num_blocks){
    oName <- omics_names[[i]] # omic names
    fNames <- names(data_blocks[[i]]) # feature names
    
    # Error catching for no names - creates default values
    if(is.null(fNames) || nchar(fNames[1])==0){
      fNames <- paste("feature",1:dim(data_blocks[[i]])[[2]],sep="_") 
    }
    
    # Checking if omics names are already at the end of feature names
    lastchars <- strsplit(fNames[[1]],split="_")
    lastchars <- lastchars[[1]][[length(lastchars[[1]])]]
    
    # If features do not have omics name at end of name, add it
    if(!tolower(lastchars) == oName & !lastchars ==oName){
      new_names <- paste(fNames,oName,sep = "_")
      names(data_blocks[[i]]) <- new_names
    }
  }
  
  
  # Pre-processing data
  if(tolower(preprocMethod) == 'colprofile'){
    message("Performing centered column profile pre-processing...")
    preproc_results <- lapply(data_blocks,CCpreproc)
    data_blocks <- lapply(preproc_results, `[[`, 1) # normalized data matrices
    block_vars <- lapply(preproc_results, `[[`, 2) # list of block variances
    message("Pre-processing completed.")
  }else{
    # **PLACEHOLDER** == should be replaced with appropriate variance calc
    data_blocks <- lapply(data_blocks, as.matrix) # converting input data to matrix form
    block_vars <- NULL
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
  eigvals <- list(nipals_result$eigval)
  
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
                      block_scores, block_loadings, eigvals, tolower(preprocMethod),
                      block_vars)
  names(results_list) <- c('global_scores','global_loadings','block_score_weights',
                           'block_scores','block_loadings', 'eigvals','preprocMethod'
                           , "block_variances")
  results_list$metadata <- metadata
  
  
  # Plotting results
  if(tolower(plots) == 'all'){
    par(mfrow = c(1,2))
    projection_plot(results_list,'projection',
               legend_loc = "bottomleft",
               color_col = color_col) # first two orders of scores
    global_scores_eigenvalues_plot(results_list) # global score eigenvalues
    par(mfrow = c(1,1))
    
  }else if (tolower(plots) == 'global'){
    par(mfrow = c(1,2))
    projection_plot(results_list,'projection_global', color_col = color_col) # first two global scores
    global_scores_eigenvalues_plot(results_list) # global score eigenvalues
    par(mfrow = c(1,1))
  }else if (tolower(plots) == 'none'){
    
  }else{
    message("No known plotting options specified - skipping plots.")
  }
  
return(results_list)
}
