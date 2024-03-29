
#' Main NIPALS computation loop
#'
#' @description Applies the full adjusted NIPALS algorithm to generate block
#' and global scores/loadings with the desired deflation method.
#'
#' @details Follows the NIPALS algorithm as described by Hanafi et. al. (2010).
#' For each order of scores/loadings, the vectors are computed via the
#' `nipals_iter` function, then used to deflate the data matrix according to
#' the desired deflation method.
#' This process is repeated up to the desired maximum order of scores/loadings.
#'
#' @param data_blocks_mae a MultiAssayExperiment class
#' object (with sample metadata as a dataframe in the colData attribute).
#' @param col_preproc_method an option for the desired column-level data
#' pre-processing, either:
#' \itemize{
#' \item `colprofile` applies column-centering, row and column weighting by
#' contribution to variance.
#' \item `standardized` centers each column and divides by its standard
#' deviation.
#' \item `centered_only` ONLY centers data
#' }
#' @param block_preproc_method an option for the desired block-level data
#' pre-processing, either:\itemize{
#' \item `unit_var` FOR CENTERED MATRICES ONLY - divides each block by the
#' square root of its variance
#' \item `num_cols` divides each block by the number of variables in the block.
#' \item `largest_sv` divides each block by its largest singular value.
#' \item `none` performs no preprocessing
#' }
#' @param num_PCs the maximum order of scores/loadings
#' @param tol a number for the tolerance on the stopping criterion for NIPALS
#' @param max_iter a number for the maximum number of times NIPALS should
#' iterate
#' @param color_col Optional argument with the column name of the `metadata`
#' data frame used to define plotting colors
#' @param deflationMethod an option for the desired deflation method, either:
#' \itemize{
#' \item `block` deflation via block loadings (for MCIA, default)
#' \item `global` deflation via global scores (for CPCA)
#' }
#' @param plots an option to display various plots of results: \itemize{
#' \item `all` displays plots of block scores, global scores, and eigenvalue
#' scree plot
#' \item `global` displays only global score projections and eigenvalue scree
#' plot
#' \item `none` does not display plots
#' }
#' @return a `nipalsResult` object with the following fields: \itemize{
#' \item `global_scores` a matrix containing global scores as columns
#' (NOT normalized to unit variance)
#' \item `global_loadings` a matrix containing global loadings as columns
#' \item `global_score_weights` a matrix of weights to express global scores as
#' a combination of block scores. Has dimensions "num_Blocks" by "num_PCs"
#' \item `eigvals` a matrix containing the eigenvalue for each computed global
#' score.
#' \item `block scores` a list of matrices, each contains the scores for one
#' block
#' \item `block loadings` a list of matrices, each contains the loadings for
#' one block (w/ unit length)
#' \item `block score weights` a matrix containing weights for each block score
#' of each order used to construct the global scores.
#' \item `col_preproc_method` the column preprocessing method used on the data.
#' \item `block_variances` a list of variances of each block AFTER
#' NORMALIZATION OPTION APPLIED
#' \item `metadata` the metadata dataframe supplied with the `metadata`
#' argument. Note: overrides metadata present in any MAE class object.
#'}
#' @importFrom graphics par
#' @importFrom methods new
#' @importFrom MultiAssayExperiment colData
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @examples
#'    data(NCI60)
#'    data_blocks_mae <- simple_mae(data_blocks,row_format="sample",
#'                                  colData=metadata_NCI60)
#'    NIPALS_results <- nipals_multiblock(data_blocks_mae, num_PCs = 10,
#'                                        tol = 1e-12, max_iter = 1000,
#'                                        col_preproc_method = "colprofile",
#'                                        deflationMethod = "block")
#'    MCIA_results <- nipals_multiblock(data_blocks_mae, num_PCs = 4)
#'    CPCA_results <- nipals_multiblock(data_blocks_mae, num_PCs = 4,
#'    deflationMethod = 'global')
#'
#' @export
nipals_multiblock <- function(data_blocks_mae,
                              col_preproc_method = "colprofile",
                              block_preproc_method = "unit_var",
                              num_PCs = 10, tol = 1e-9,
                              max_iter = 1000,color_col = NULL,
                              deflationMethod = "block", plots = "all") {
    # Check for input type MAE or list
    if (is(data_blocks_mae, "MultiAssayExperiment")) {
      data_blocks <- extract_from_mae(data_blocks_mae)
      metadata <- data.frame(MultiAssayExperiment::colData(data_blocks_mae))
      if (length(metadata) == 0) {
        metadata <- NULL
      }
    } else {
        stop("Unknown input data format - please use MultiAssayExperiment")
    }

    num_blocks <- length(data_blocks)
    omics_names <- names(data_blocks)

    # Check for omics names and assign generic name if null
    if (is.null(omics_names)) {
        omics_names <- paste("omic", seq(1, length(data_blocks)), sep = "")
        names(data_blocks) <- omics_names
    }

    # Check for feature names and assign generic name if null
    for (i in seq_along(data_blocks)) {
        fNames <- names(data_blocks[[i]]) # feature names

        # Error catching for no names - creates default values
        if (is.null(fNames) || nchar(fNames[1]) == 0) {
            fNames <- paste("feature", seq(1, dim(data_blocks[[i]])[[2]]),
                            sep = "_")
        }
    }

    # Pre-processing data
    message("Performing column-level pre-processing...")
    data_blocks <- lapply(data_blocks, col_preproc, col_preproc_method)
    message("Column pre-processing completed.")

    # Block-level pre-processing
    message("Performing block-level preprocessing...")
    data_blocks <- lapply(data_blocks, block_preproc, block_preproc_method)
    if (tolower(block_preproc_method) == "unit_var") {
        block_vars <- rep(list(1), length(data_blocks))
        names(block_vars) <- names(data_blocks)

    } else {
        block_vars <- get_tv(data_blocks)
    }

    message("Block pre-processing completed.")

    # First NIPALS run
    message("Computing order ", 1, " scores")
    nipals_result <- nipals_iter(data_blocks, tol)

    # Saving result
    # matrix containing global scores as columns
    global_scores <- nipals_result$global_scores
    # matrix containing global loadings as columns
    global_loadings <- as.matrix(nipals_result$global_loadings)
    # matrix containing block score weights as columns
    block_score_weights <- nipals_result$block_score_weights

    block_scores <- list() # list containing matrices of block scores
    block_loadings <- list() # list containing matrices of block loadings
    for (i in seq_along(data_blocks)) {
        block_scores[[i]] <- nipals_result$block_scores[, i]
        block_loadings[[i]] <- nipals_result$block_loadings[[i]]
    }

    # Computing block eigenvalue
    eigvals <- list(nipals_result$eigval)

    if (num_PCs > 1) {
        # generate scores/loadings up to number of PCs
        for (i in seq(2, num_PCs)) {
            message("Computing order ", i, " scores")
            # Deflate blocks
            if (tolower(deflationMethod) == "block") {
                data_blocks <- mapply(deflate_block_bl,
                                      data_blocks,
                                      nipals_result$block_loadings,
                                      SIMPLIFY = FALSE)

            } else if (tolower(deflationMethod) == "global") {
                data_blocks <- lapply(data_blocks,
                                      deflate_block_gs,
                                      gs = nipals_result$global_scores)

            } else {
                stop("Unknown option for deflation step -
                     use 'block' or 'global'")
            }

            # Run another NIPALS iteration
            nipals_result <- nipals_iter(data_blocks, tol)

            # Save results
            global_scores <- cbind(global_scores, nipals_result$global_scores)
            global_loadings <- cbind(global_loadings,
                                     nipals_result$global_loadings)
            block_score_weights <- cbind(block_score_weights,
                                         nipals_result$block_score_weights)
            eigvals <- cbind(eigvals, nipals_result$eigval)

            for (j in seq(1, num_blocks)) {
                block_scores[[j]] <- cbind(block_scores[[j]],
                                           nipals_result$block_scores[, j])
                block_loadings[[j]] <- cbind(block_loadings[[j]],
                                             nipals_result$block_loadings[[j]])
            }
        }
    }

    # Formatting results
    # adding block names in outputs
    names(block_scores) <- names(data_blocks)
    names(block_loadings) <- names(data_blocks)
    names(eigvals) <- paste("gs", seq(1, num_PCs), sep = "")

    # adding row (sample) names in outputs
    rownames(global_scores) <- rownames(data_blocks[[1]])

    for (j in seq(1, num_blocks)) {
        rownames(block_scores[[j]]) <- rownames(data_blocks[[1]])
    }

    # fixing no metadata with empty data frame
    if (is.null(metadata)) {
        metadata <- data.frame()
    }

    # creating S4 class object with nipals outputs
    mcia_out <- new("NipalsResult",
                    global_scores = global_scores,
                    global_loadings = global_loadings,
                    block_score_weights = block_score_weights,
                    block_scores = block_scores,
                    block_loadings = block_loadings,
                    eigvals = eigvals,
                    col_preproc_method = tolower(col_preproc_method),
                    block_preproc_method = tolower(block_preproc_method),
                    block_variances = block_vars,
                    metadata = metadata)

    # Plotting results
    if (tolower(plots) == "all") {
        par(mfrow = c(1, 2))
        projection_plot(mcia_out, "all", legend_loc = "bottomleft",
                        color_col = color_col) # first two orders of scores
        global_scores_eigenvalues_plot(mcia_out) # global score eigenvalues
        par(mfrow = c(1, 1))

    } else if (tolower(plots) == "global") {
        par(mfrow = c(1, 2))
        # first two global scores
        projection_plot(mcia_out, "global", color_col = color_col)
        global_scores_eigenvalues_plot(mcia_out) # global score eigenvalues
        par(mfrow = c(1, 1))

    } else if (tolower(plots) == "none") {

    } else {
        message("No known plotting options specified - skipping plots.")
    }

    return(mcia_out)
}
