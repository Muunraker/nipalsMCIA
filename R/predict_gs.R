
#' Prediction of new global scores based on block loadings and weights
#'
#' @description Uses previously-computed block scores and weights to compute
#' a global score for new data.
#'
#' @details Projects the new observations onto each block loadings vector, then
#' weights the projection according to the corresponding block weights.
#'
#' @param mcia_results an mcia object output by nipals_multiblock()
#'         containing block scores, weights, and pre-processing identifier.
#' @param df a list of data matrices to make predictions from, where each entry
#'         corresponds to one omics type in "sample" x " features" format.
#'         Feature and omic order must match `bl`. Pre-processing should also
#'         match the data used to generate `bl` and `bw`.
#' @return a matrix of predicted global scores, in form
#' @examples
#'    data(NCI60)
#'    data_blocks_mae <- simple_mae(data_blocks,row_format="sample",
#'                                  colData=metadata_NCI60)
#'    mcia_results <- nipals_multiblock(data_blocks_mae, num_PCs = 2)
#'    new_data <- data_blocks # should update with a truly new dataset
#'    preds <- predict_gs(mcia_results, new_data)
#'
#' @export
#'
predict_gs <- function(mcia_results, df) {
    bl <- mcia_results@block_loadings
    bw <- mcia_results@block_score_weights
    col_preproc_method <- mcia_results@col_preproc_method
    block_preproc_method <- mcia_results@block_preproc_method

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

    # if (tolower(preproc_method) == "colprofile") {
    #    message("Centered column profile pre-processing detected.")
    #    message("Performing pre-processing on data")
    #    df <- lapply(df, CCpreproc)
    #    message("Pre-processing completed.")
    # } else {
    #    message("No pre-processing detected.")
    #    df <- lapply(df, as.matrix)
    # }

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
