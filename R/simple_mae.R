#' Create an MAE object from a list of data matrices and column data
#'
#' @description Create an MAE object from a set of data matrices and column data.
#' 
#' @details Requires that sample names match across experiments and are identical
#' to primary names, will only convert data matrices to SummarizedExperiment class.
#' If the data is more complex, please follow the guidelines 
#' for creating an MAE object outlined in `help(MultiAssayExperiment)`
#'
#' @param matrix_list named list of data matrices
#' @param row_format for lists of data frames, indicates whether rows of
#' datasets denote `feature` (default) or `sample`.
#' @param colData_input a data frame containing sample metadata; sample names
#' in the rownames should correspond to samples names in `matrix_list`
#' @return List of harmonized data matrices for input into nipals_multiblock()
#' @examples
#' data(NCI60)
#' data_blocks_mae <- simple_mae(data_blocks,row_format='sample', 
#'                               colData=metadata_NCI60)
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment 
#' @export


simple_mae <- function(matrix_list,row_format = "feature", colData_input=NULL){
  
  # The MAE class takes in data matrices in feature x sample orientation
  if (row_format=="sample"){matrix_list<-lapply(matrix_list, t)}
  else if (row_format!="feature"){
    stop("Choose from sample or feature orientation")
  }
  
  if(is.null(colData_input)){
  mae_object<-MultiAssayExperiment::MultiAssayExperiment(
      lapply(matrix_list, function(x) SummarizedExperiment(as.matrix(x))))}
  else {mae_object<-MultiAssayExperiment::MultiAssayExperiment(
    lapply(matrix_list, function(x) SummarizedExperiment(as.matrix(x))),
    colData = colData_input)}

return(mae_object)
}
