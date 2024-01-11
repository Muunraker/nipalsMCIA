#' An S4 class to contain results computed with `nipals_multiblock()`
#' @slot global_scores A matrix containing global scores as columns.
#' @slot global_loadings A matrix containing global loadings as columns.
#' @slot block_score_weights A matrix containing block weights as columns.
#' @slot block_scores A list of matrices. Each matrix contains the scores
#' as columns for a given block.
#' @slot block_loadings A list of matrices. Each matrix contains the
#' loadings as columns for a given block.
#' @slot eigvals A list of singular values of the data matrix at each
#' deflation step.
#' @slot col_preproc_method character for the column-level preprocessing
#' method used. See `col_preproc()`.
#' @slot block_preproc_method character for the block-level
#' preprocessing method used. See `block_preproc()`.
#' @slot block_variances A list of variances for each block.
#' @slot metadata A data frame of metadata originally passed into
#' `nipals_multiblock()`.
#'
#' @returns A NipalsResult object.
#' @exportClass NipalsResult

NipalsResult <- setClass("NipalsResult",
  representation(
    global_scores = "matrix",
    global_loadings = "matrix",
    block_score_weights = "matrix",
    block_scores = "list",
    block_loadings = "list",
    eigvals = "matrix",
    col_preproc_method = "character",
    block_preproc_method = "character",
    block_variances = "list",
    metadata = "data.frame")
)

# defining show method for NipalsResult
setMethod("show", "NipalsResult",
  function(object){
    cat(is(object)[[1]], " Object with properties: \n",
      "> Dataset dimensions:   ", nrow(object@global_scores), " x ", nrow(object@global_loadings), "\n",
      "> Number of blocks:  ", length(object@block_variances), "\n",
      "> Order of scores:  ", ncol(object@global_scores), "\n",
      "> Column preprocessing:  ", object@col_preproc_method, "\n",
      "> Block preprocessing:  ", object@block_preproc_method, "\n",
      "> Block names and sizes:   \n",
      sep = ""
    )
    blockdims <- lapply(object@block_loadings, nrow)
    print(unlist(blockdims))
  }
)
