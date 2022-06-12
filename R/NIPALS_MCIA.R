

#' Omicade4 Initialization
#'
#' @description Applies the variable and block level pre-processing as the 
#' "Omicade4" package (Meng et al. 2014)
#' @details Performs the following steps on a given data frame: \itemize{
#' \item Offsets data to make whole matrix non-negative
#' \item Divides each column by its sum 
#' \item Subtracts (row sum/total sum) from each row
#' \item Multiplies each column by sqrt(column sum/total sum)
#' \item Divides the whole data frame by its total variance (the sqrt of the sum of singular values)
#' }
#' @param df the data frame to apply pre-processing to, in "sample" x "variable" format 
#' @return the processed data frame
#' @export
omicadeInitialization <- function(df){
  temp_df <- df
  
  # Making data non-negative
  minVal <- min(temp_df)
  if(minVal < 0 ){
    offset <- floor(minVal)
    temp_df <- temp_df+offset
  }

  # Generating centered column profiles:
  totsum <- sum(temp_df)
  colsums <- colSums(temp_df)
  row_contribs <- rowSums(temp_df)/totsum
  
  ## Dividing by column sums
  nz_cols <- which(colsums != 0) # excluding zero columns to avoid NaNs 
  temp_df[,nz_cols] <- as.data.frame(t( t(temp_df[,nz_cols])/colsums[nz_cols]))
  
  ## Subtracting row contributions
  temp_df <- temp_df - row_contribs

  # Applying feature weighting ("multiplication by feature metrics")
  temp_df <- as.data.frame(t( t(temp_df)*sqrt(colsums/totsum)))
  
  # Applying block weights (blocks have unit variance via division by sum of eigenvalues)
  temp_df <- temp_df*(1/sqrt(sum(svd(temp_df)[[1]]^2)))
  
  return(temp_df)
}

#' NIPALS Iteration
#'
#' @description Applies one iteration stage/loop of the NIPALS algorithm. 
#' 
#' @details Follows the NIPALS algorithm as described by Hanafi et. al. (2010).
#' Starts with a random vector in sample space and repeatedly projects it onto
#' the variable vectors and block scores to generate block and global 
#' loadings/scores/weights. The loop stops when either the stopping criterion is 
#' low enough, or the maximum number of iterations is reached. Intended as a 
#' utility function for `nipals_multiblock` to be used between deflation steps.
#' 
#' @param ds a list of data frames, each in "sample" x "variable" format  
#' @param tol a number for the tolerance on the stopping criterion for NIPALS
#' @param maxIter a number for the maximum number of times NIPALS should iterate
#' @return a list containing the global/block scores, loadings and weights for a given order
#' @examples 
#' nipals_results <- NIPALS_iter(data_list, tol = 1e-7, maxIter = 1000)
#' 
#' @export
NIPALS_iter <- function(ds, tol=1e-10, maxIter=1000){
  
  # Main iteration loop
  stopCrit <- 2*tol  
  covSquared_old <- 0
  iter <- 0
  gs <- pracma::rand(nrow(ds[[1]]),1) # begin with random global score vector
  
  while(stopCrit > tol && iter <= maxIter){
    
    # Computing block loadings
    bl_list <- lapply(ds, function(df,q){ 
      bl_k <- crossprod(as.matrix(df), q) 
      bl_k <- bl_k/norm(bl_k, type="2") 
      return(bl_k)
    }, q=gs)
    
    # Computing block scores
    bs_list <- mapply(function(df,bl_k){
      bs_k <- as.matrix(df) %*% bl_k
      return(bs_k)
    },ds, bl_list)
    
    # Computing global weights
    gw <- crossprod(bs_list,gs)
    gw <- gw/norm(gw, type="2")
    gs <- bs_list %*% gw
    
    # Computing stopping criteria
    covList <- sapply(as.data.frame(bs_list), function(bs, gs){
      gs_norm <- gs/sqrt(drop(var(gs)))
      return(drop(cov(bs,gs_norm))^2)
    }, gs = gs)
    
    stopCrit <- abs(sum(covList) - covSquared_old)
    covSquared_old <-sum(covList)
    
    iter <- iter +1 
    
    print(iter)
    print(stopCrit)
    
    
  }
  if(iter > maxIter){
    warning('NIPALS iteration did not converge')
  }
  
  # Computing global loadings at final iteration
  gl <- bl_list[[1]]*gw[1]
  nblocks <- length(ds)
  for(i in 2:nblocks){
    gl <- c(gl, bl_list[[i]]*gw[i])
  }
  
  # Returning results
  retlist <-list(gs, gl, gw, bs_list, bl_list )
  names(retlist) <- c('global_scores','global_loadings','block_score_weights','block_scores','block_loadings')
  return(retlist)
}


#' Deflation via block loadings
#'
#' @description Removes data from a data frame in the direction of a given block 
#' loadings vector. 
#' 
#' @details Subtracts the component of each row in the direction of a given
#' block loadings vector to yield a `deflated' data matrix.
#' 
#' @param df a data frame in "sample" x "variable" format  
#' @param bl a block loadings vector in variable space
#' @return the deflated data frame
#' @examples 
#' deflated_data <- deflate_block_bl(data_frame,block_loading)
#' 
#' @export
deflate_block_bl <- function(df,bl){
  df <- df - tcrossprod(as.matrix(df) %*% bl, bl)
  return(df)
}

#' Deflation via global scores
#'
#' @description Removes data from a data frame in the direction of a given global 
#' scores vector. 
#' 
#' @details Subtracts the component of each column in the direction of a given
#' global scores vector to yield a `deflated' data matrix.
#' 
#' @param df a data frame in "sample" x "variable" format  
#' @param gs a global scores vector in sample space
#' @return the deflated data frame
#' @examples 
#' deflated_data <- deflate_block_gs(data_frame,global_score)
#' 
#' @export
deflate_block_gs <- function(df,gs){
  normed_gs <- gs/norm(gs, type="2")
  df <- df - tcrossprod(normed_gs) %*% as.matrix(df)
  return(df)
}


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
#' @param num_PCs the maximum order of scores/loadings
#' @param tol a number for the tolerance on the stopping criterion for NIPALS
#' @param maxIter a number for the maximum number of times NIPALS should iterate
#' @param deflationMethod an option for the desired deflation method, either: \itemize{
#' \item `block` deflation via block loadings (for MCIA, default)
#' \item `global` deflation via global scores (for CPCA)
#' }
#' @return a list containing the following: \itemize{
#' \item `global_scores` a matrix containing global scores as columns
#' \item `global_loadings` a matrix containing global loadings as columns
#' \item `global_score_weights` a matrix of weights to express global scores as
#' a combination of block scores. Has dimensions "num_Blocks" by "num_PCs"
#' \item `block scores` a list of matrices, each contains the scores for one block
#' \item `block loadings` a list of matrices, each contains the loadings for one block
#' }
#' @examples 
#'  NIPALS_results <- nipals_multiblock(df_list, num_PCs = 2, tol = 1e-7, maxIter = 1000, deflationMethod = 'block')
#'  MCIA_result <- nipals_multiblock(df_list, num_PCs = 2)
#'  CPCA_result <- nipals_multiblock(df_list, num_PCs = 4,deflationMethod = 'global')
#' 
#' @export
nipals_multiblock <- function(data_blocks, num_PCs, tol=1e-7, max_iter = 1000, deflationMethod = 'block'){
  num_blocks <- length(data_blocks)
  
  # First NIPALS run
  nipals_result <- NIPALS_iter(data_blocks, tol)
  
  # Setting up variables
  global_scores <- nipals_result$global_scores # matrix containing global scores as columns
  global_loadings <- as.matrix(nipals_result$global_loadings) # matrix containing global loadings as columns
  block_score_weights <- nipals_result$block_score_weights # matrix containing block score weights as columns
  
  block_scores <- list() # list containing matrices of block scores
  block_loadings <- list() # list containing matrices of block loadings
  for(i in 1:num_blocks){
    block_scores[[i]] <- nipals_result$block_scores[,i]
    block_loadings[[i]] <- nipals_result$block_loadings[[i]]
  }
                       

  # generate scores/loadings up to number of PCs
  for(i in 2:num_PCs){
    # Deflate blocks
    if(tolower(deflationMethod) == 'block'){
      data_blocks <- mapply(deflate_block_bl, data_blocks, nipals_result$block_loadings)
    } else if(tolower(deflationMethod) == 'global'){
      data_blocks <- lapply(data_blocks, deflate_block_gs, gs=nipals_result$global_scores)
    }else{
      stop("Uknown option for deflation step - use 'block' or 'global'")
    }
    
    nipals_result <- NIPALS_iter(data_blocks, tol)
    
    global_scores <- cbind(global_scores, nipals_result$global_scores)
    global_loadings <- cbind(global_loadings, nipals_result$global_loadings)
    block_score_weights <- cbind(block_score_weights, nipals_result$block_score_weights)
    
    for(j in 1:num_blocks){
      block_scores[[j]] <- cbind(block_scores[[j]], nipals_result$block_scores[,j])
      block_loadings[[j]] <- cbind(block_loadings[[j]], nipals_result$block_loadings[[j]])
    }
  }
  names(block_scores) <- names(data_blocks)
  names(block_loadings) <- names(data_blocks)
  retlist <-list(global_scores, global_loadings, block_score_weights, block_scores, block_loadings )
  names(retlist) <- c('global_scores','global_loadings','block_score_weights','block_scores','block_loadings')
  return(retlist)
}
