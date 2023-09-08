#' Ranked global loadings dataframe
#'
#' @description Creates a dataframe with ranked loadings for a given factor
#'
#' @param mcia_out object returned from nipals_multiblock() function
#' @param omic choose an omic to rank, or choose 'all' for all,
#' ((omic = "all", omic = "miRNA", etc.))
#' @param absolute whether to rank by absolute value
#' @param descending whether to rank in descending or ascending order
#' @param factor choose a factor (numeric value from 1 to number of factors
#'               in mcia_out)
#' @return ranked dataframe
#' @examples
#' data(NCI60)
#' data_blocks_mae <- simple_mae(data_blocks,row_format="sample",
#'                               colData=metadata_NCI60)
#' mcia_results <- nipals_multiblock(data_blocks_mae, num_PCs = 10,
#'                                   plots = "none", tol = 1e-12)
#' all_pos_1 <- ord_loadings(mcia_out = mcia_results, omic = "all",
#'                           absolute = FALSE, descending = TRUE, factor = 1)
#' @export

ord_loadings <- function(mcia_out, omic = "all", absolute = FALSE,
                         descending = TRUE, factor = 1) {
    # list of omics plus 'all'
    omics_names <- names(mcia_out@block_loadings)
    omics_names <- c(omics_names, "all")

    # Return error if omic not chosen from omics_names list
    if (!(omic %in% omics_names)) {
        stop("Choose an appropriate omic")
    }

    # Get global loadings and list of omics
    gl <- mcia_out@global_loadings
    # omic_type <- gsub("^.*_", "", rownames(gl))

    omic_dims <- vapply(mcia_out@block_loadings, dim, numeric(2))[1, ]
    omic_type <- c()
    omics_labels <- names(mcia_out@block_loadings)
    for (i in seq_along(omics_labels)) {
        omic_label <- omics_labels[i]
        length_omic <- omic_dims[i]
        omic_type <- c(omic_type, rep(omic_label, each = length_omic))
    }

    # Filter on factor, and add omic type
    gl_f <- data.frame(gl[, factor])
    gl_f$omic <- omic_type
    gl_f$omic <- as.factor(gl_f$omic)
    colnames(gl_f) <- c("loading", "omic")

    # Rank features by parameter settings
    if (absolute == FALSE) {
        if (descending == TRUE) {
            gl_f_ord <- gl_f[order(gl_f$loading, decreasing = TRUE), ]
        } else {
            gl_f_ord <- gl_f[order(gl_f$loading, decreasing = FALSE), ]
        }

    } else {
        gl_f_abs <- gl_f
        gl_f_abs$abs <- abs(gl_f_abs$loading)
        gl_f_ord <- gl_f_abs[order(gl_f_abs$abs, decreasing = TRUE), ]
    }

    if (omic != "all") {
        gl_f_ord <- gl_f_ord[gl_f_ord$omic == omic, ]
    }

    gl_f_ord$omic_name <- rownames(gl_f_ord)
    gl_f_ord$factor <- factor

    return(gl_f_ord)
}
