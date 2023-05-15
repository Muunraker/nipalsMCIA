#' Computes the total variance of a multi-omics dataset
#'
#' @description Computes the total variances of all data blocks in a
#' multi-omics dataset, intended for datasets that do not use `CCpreproc`
#' @param ds a list of multi-omics dataframes/matrices in "sample x variable"
#' format
#' @return the total variance of the dataset (i.e. sum of block variances)
#' @examples
#' data(NCI60)
#' tot_var <- get_TV(data_blocks)
#' @export
get_tv <- function(ds) {
    ds <- lapply(ds, as.matrix)

    # function to calculate variances of each block
    get_bv <- function(X) {
        # comment out for centered matrices
        X_cent <- sweep(X, 2, colMeans(X), "-")

        n_samp <- dim(X_cent)[1]
        X_cent_norm <- X_cent

        bv <- (norm(X_cent_norm, type = "F")^2) / max(n_samp - 1, 1)
        return(bv)
    }

    bv_list <- lapply(ds, get_bv)
    # tv <- sum(unlist(bv_list))

    return(bv_list)
}
