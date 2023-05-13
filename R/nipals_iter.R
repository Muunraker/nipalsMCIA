#' NIPALS Iteration
#'
#' @description Applies one iteration stage/loop of the NIPALS algorithm.
#'
#' @details Follows the NIPALS algorithm as described by Hanafi et. al. (2010).
#' Starts with a random vector in sample space and repeatedly projects it onto
#' the variable vectors and block scores to generate block and global
#' loadings/scores/weights. The loop stops when either the stopping criterion
#' is low enough, or the maximum number of iterations is reached. Intended as a
#' utility function for `nipals_multiblock` to be used between deflation steps.
#'
#' @param ds a list of data matrices, each in "sample" x "variable" format
#' @param tol a number for the tolerance on the stopping criterion for NIPALS
#' @param maxIter a number for the maximum number of times NIPALS should iterate
#' @return a list containing the global/block scores, loadings and weights for
#' a given order
#' @examples
#' data(NCI60)
#' data_blocks <- lapply(data_blocks, as.matrix)
#' nipals_results <- NIPALS_iter(data_blocks, tol = 1e-7, maxIter = 1000)
#' @importFrom pracma rand
#' @importFrom RSpectra svds
#' @importFrom stats cov var
#' @export
nipals_iter <- function(ds, tol = 1e-12, maxIter = 1000) {
    # Main iteration loop
    stop_crit <- 2 * tol
    cov_squared_old <- 0
    iter <- 0
    gs <- pracma::rand(nrow(ds[[1]]), 1) # begin with random global score vector

    while (stop_crit > tol && iter <= maxIter) {
        # Computing block loadings
        bl_list <- lapply(ds, function(df, q) {
            bl_k <- crossprod(df, q)
            bl_k <- bl_k / norm(bl_k, type = "2")
            return(bl_k)
        }, q = gs)

        # Computing block scores
        bs_list <- mapply(function(df, bl_k) {
            bs_k <- df %*% bl_k
            return(bs_k)
        }, ds, bl_list)

        # Computing global weights
        gw <- crossprod(bs_list, gs)
        gw <- gw / norm(gw, type = "2")
        gs <- bs_list %*% gw

        # Computing stopping criteria
        cov_list <- vapply(as.data.frame(bs_list), function(bs, gs) {
            gs_norm <- gs / sqrt(drop(var(gs)))
            return(drop(cov(bs, gs_norm))^2)
        }, gs = gs, FUN.VALUE = numeric(1))

        stop_crit <- abs(sum(cov_list) - cov_squared_old)
        cov_squared_old <- sum(cov_list)

        iter <- iter + 1
    }

    if (iter > maxIter) {
        warning("NIPALS iteration did not converge")
    }

    # Computing eigenvalue associated with the global score
    global_matrix <- do.call(cbind, ds)
    svdres <- RSpectra::svds(global_matrix, 1)
    eigval <- (svdres$d)^2 / (max(dim(global_matrix)[1] - 1, 1))

    # Computing global loadings at final iteration
    gl <- bl_list[[1]] * gw[1]
    names(gl) <- rownames(bl_list[[1]])
    nblocks <- length(ds)
    for (i in seq(2, nblocks)) {
        temp_names <- c(names(gl), rownames(bl_list[[i]]))
        gl <- c(gl, bl_list[[i]] * gw[i])
        names(gl) <- temp_names
    }

    # Returning results
    retlist <- list(gs, gl, gw, bs_list, bl_list, eigval)
    names(retlist) <- c("global_scores", "global_loadings",
                        "block_score_weights", "block_scores",
                        "block_loadings", "eigval")
    return(retlist)
}
