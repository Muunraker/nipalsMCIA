
#' Prediction of new global scores based on block loadings and weights
#'
#' @description Uses previously-computed block scores and weights to compute 
#' a global score for new data.
#' 
#' @details Projects the new observations onto each block loadings vector, then 
#' weights the projection according to the corresponding block weights.
#' 
#' @param bl a list of matrices of block loadings, where each matrix corresponds
#' to one omics type with dimensions "features" x "number of loadings" 
#' @param bw a matrix of block weights, with dimensions "number of omics" x "number of loadings"
#' @param df a list of data matrices to make predictions from, where each entry 
#' corresponds to one omics type in "sample" x " features" format.
#' Feature and omic order must match `bl`. Pre-processing should also match the data
#' used to generate `bl` and `bw`. 
#' @return a matrix of predicted global scores, in form 
#' @examples 
#' deflated_data <- deflate_block_bl(data_frame,block_loading)
#' 
#' @export
#' 
predict_gs <- function(bl,bw, df){
  
  num_omics <- length(bl)
  if(length(df) != length(bl) | length(df) != dim(bw)[[1]]){
    stop("Mismatched number of omics between the block loadings, block weights, and new data.")
  }
  
  
  # ensuring all arguments are matrices
  bl <- lapply(bl,as.matrix)
  df <- lapply(df,as.matrix)
  
  if(dim(df[[1]])[[2]] != dim(bl[[1]])[[1]]){
    stop(paste('Error: mismatched number of features in omic ',1))
  }
  new_gs <- df[[1]] %*% bl[[1]] # block score matrix for 1st omic
  new_gs <- t(t(new_gs)*bw[1,]) # applying block weight
  
  # for each omics type, a
  if(num_omics >1){
    for( i in 2:num_omics){
      if(dim(df[[i]])[[2]] != dim(bl[[i]])[[1]]){
        stop(paste('Error: mismatched number of features in omic ',i))
      }
      new_gs_i <- df[[i]] %*% bl[[i]] # block score matrix for ith omic
      new_gs <- new_gs +t(t(new_gs_i)*bw[i,]) # applying block weight
    }
  }
  return(new_gs)
}
