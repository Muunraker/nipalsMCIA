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
#'     in mcia_out)
#' @return ranked dataframe
#' @examples
#' all_pos_1 <- ord_loadings(mcia_out = mcia_results, omic = "all",
#'                          absolute = FALSE, descending = TRUE, factor = 1)
#' @export

ord_loadings <- function(mcia_out, omic = "all", absolute = FALSE,
                         descending = TRUE, factor = 1) {
  # list of omics plus 'all'
  omics_names <- names(mcia_out$block_loadings)
  omics_names <- c(omics_names, "all")

  # Return error if omic not chosen from omics_names list
  if (!(omic %in% omics_names)) {
      stop("Error: choose appropriate omic")
  }

  # Get global loadings and list of omics
  gl <- mcia_out$global_loadings
  omic_type <- gsub("^.*_", "", rownames(gl))

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
