
#' Prediction of new global scores based on block loadings and weights
#'
#' @description Uses previously-computed block scores and weights to compute
#' a global score for new data. Only validated for MCIA results, as CPCA loadings
#' aren't compatible with un-deflated data.
#'
#' @details Projects the new observations onto each block loadings vector, then
#' weights the projection according to the corresponding block weights.
#'
#' @param mcia_results an mcia object output by nipals_multiblock()
#'         containing block scores, weights, and pre-processing identifier.
#' @param test_data an MAE object with the same block types and features as the
#'         training dataset. Feature and omic order must match `bl`.
#' @return a matrix of predicted global scores for the training data
#' @examples
#'    data(NCI60)
#'    data_blocks_mae <- simple_mae(data_blocks,row_format="sample",
#'                                  colData=metadata_NCI60)
#'    mcia_results <- nipals_multiblock(data_blocks_mae, num_PCs = 2)
#'    new_data <- data_blocks_mae # should update with a truly new dataset
#'    preds <- predict_gs(mcia_results, new_data)
#'
#' @export
#'
predict_gs <- function(mcia_results, test_data) {
    bl <- mcia_results@block_loadings
    bw <- mcia_results@block_score_weights
    col_preproc_method <- mcia_results@col_preproc_method
    block_preproc_method <- mcia_results@block_preproc_method

    # extract data from train MAE object
    df <- extract_from_mae(test_data)

    num_omics <- length(bl)
    if (length(df) != length(bl) || length(df) != dim(bw)[[1]]) {
        stop("Mismatched number of omics between the block loadings, ",
             "block weights, and new data.")
    }

    # Apply same pre-processing methods as the model
    message("Performing pre-processing on data")
    df <- lapply(df, col_preproc, col_preproc_method)
    df <- lapply(df, block_preproc, block_preproc_method)
    message("Pre-processing completed")

    # ensuring all arguments are matrices
    bl <- lapply(bl, as.matrix)

    if (dim(df[[1]])[[2]] != dim(bl[[1]])[[1]]) {
        stop("Mismatched number of features in omic ", 1)
    }

    new_gs <- df[[1]] %*% bl[[1]] # block score matrix for 1st omic
    new_gs <- t(t(new_gs) * bw[1, ]) # applying block weight

    # for each omics type
    if (num_omics > 1) {
        for (i in 2:num_omics) {
            if (dim(df[[i]])[[2]] != dim(bl[[i]])[[1]]) {
                stop("Mismatched number of features in omic ", i)
            }

            new_gs_i <- df[[i]] %*% bl[[i]] # block score matrix for ith omic
            new_gs <- new_gs + t(t(new_gs_i) * bw[i, ]) # applying block weight
        }
    }

    return(new_gs)
}
