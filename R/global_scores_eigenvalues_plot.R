#' global_scores_eigenvalues_plot
#'
#' @description Function to plot eigenvalues of scores up to num_PCs
#'
#' @details Plotting function for eigenvalues of scores up to num_PCs
#'
#' @param mcia_results MCIA results object returned from `nipals_multiblock`
#'
#' @examples
#' data(NCI60)
#' mcia_results <- nipals_multiblock(data_blocks, metadata = metadata_NCI60,
#'                               num_PCs = 10, plots = "none", tol=1e-12)
#' global_scores_eigenvalues_plot(mcia_results)
#' @return Displays the contribution plot using eigenvalues
#' @importFrom graphics barplot
#' @export
global_scores_eigenvalues_plot <- function(mcia_results) {

  # Getting total variance if supplied
  if (is.list(mcia_results$block_variances)) {
    ylabel <- "Prop. Total Variance"
    totvar <- sum(unlist(mcia_results$block_variances))
  } else {
    ylabel <- "Eigenvalue"
    totvar <- 1
  }

  # extract eigenvalues
  barploteigs <- (unlist(mcia_results$eigvals)) / totvar
  names(barploteigs) <- seq(1, length(mcia_results$eigvals))

  # generate barplot
  barplot(barploteigs, xlab = "Factors", ylab = ylabel, cex.names = 1,
          main = "Global Eigenvalues ")
}
