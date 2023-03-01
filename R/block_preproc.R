#' Block-level preprocessing
#'
#' @description A function that normalizes an input dataset (data block)
#' according to a variety of options.
#' Intended to be used after column/row-level normalization.
#' @param df dataset to preprocess (must be in data matrix form)
#' @param block_preproc_method method which is used to normalize blocks, with
#' options: \itemize{
#' \item `unit_var` FOR CENTERED MATRICES ONLY - divides each block by the
#' square root of its variance
#' \item `num_cols` divides each block by the number of variables in the block.
#' \item `largest_sv` divides each block by its largest singular value.
#' \item `none` performs no preprocessing
#' }
#' @return the preprocessed dataset
#' @examples
#' df <- matrix(rbinom(15, 1, prob = 0.3), ncol = 3)
#' preprocessed_dataframe <- block_preproc(df,"unit_var")
#' @importFrom RSpectra svds
#' @export

block_preproc <- function(df, block_preproc_method) {
    temp_df <- df
    block_preproc_method <- tolower(block_preproc_method)

    # Weighting block to unit variance (block variance computed via Fro norm)
    if (block_preproc_method == "unit_var") {
        block_var <- norm(temp_df, type = "F")^2 / (max(1, nrow(temp_df) - 1))

        if (block_var != 0) {
            temp_df <- temp_df * (1 / sqrt(block_var))
        } else {
            warning("Data block has zero variance.")
        }

        return(temp_df)
    } else if (block_preproc_method == "num_cols") {
        temp_df <- temp_df * (1 / ncol(temp_df))
        return(temp_df)

    } else if (block_preproc_method == "largest_sv") {
        svdres <- RSpectra::svds(temp_df, 1)
        eigval <- svdres$d # largest singular value of the data block

        if (eigval != 0) {
            temp_df <- temp_df * (1 / eigval)
        } else {
            warning("Data block has zero largest eigenvalue")
        }

        return(temp_df)

    } else if (block_preproc_method == "none") {
        return(temp_df)

    } else {
        stop("block preprocessing method not recognized",
             "- pick from available options")
    }
}
