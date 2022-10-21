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
#' @export
global_scores_eigenvalues_plot <- function(mcia_result){
    
  # extract eigenvalues
  barploteigs <- unlist(mcia_result$eigvals)^2
  names(barploteigs) <- 1:length(mcia_result$eigvals)
  
  # generate barplot 
  barplot(barploteigs, xlab="Factors", cex.names = 1, 
          main = "Global Factor Score Eigenvalues ")
}