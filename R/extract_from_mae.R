#' Extract a list of harmonized data matrices from an MAE object
#'
#' @description Extract a list of harmonized data matrices for
#' input into nipals_multiblock() from an MAE object
#'
#' @param MAE_object an MAE object containing experiment data for extraction
#' colData field optional
#' experiments should either be SummarizedExperiment, SingleCellExperiment, or
#' RangedSummarizedExperiment classes
#' @param subset_data \itemize{
#' \item `all` use all experiments in MAE object
#' \item `c(omic1,omic2,...)` list of omics from names(MAE_object)
#' }
#' @param harmonize A boolean whether samples should be checked for duplicates \itemize{
#' \item `TRUE` (default) merges duplicate samples via the `MultiAssayExperiment::mergeReplicates` function
#' \item `FALSE` skips sample duplicate check - USE THIS FOR LARGE-SAMPLE DATASETS.
#' }
#' @return List of harmonized data matrices for input into `nipals_multiblock()`
#' @examples
#' data(NCI60)
#' data_blocks_mae <- simple_mae(data_blocks, row_format = "sample",
#'                               colData = metadata_NCI60)
#' NCI60_input = extract_from_mae(data_blocks_mae, subset = "all")
#' @importFrom MultiAssayExperiment colData assays experiments
#' @importFrom MultiAssayExperiment mergeReplicates intersectColumns
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export


extract_from_mae <- function(MAE_object, subset_data = "all", harmonize = TRUE) {
  if (subset_data != "all") {
    if (intersect(subset_data, names(MAE_object)) < 2) {
      stop("Please provide appropriate subset_data list")
    }
    else {
      MAE_object <- MAE_object[,,subset]
    }
  }

  n_omics <- length(MAE_object)
  if (n_omics == 1) {
    stop("The nipalsMCIA algorithm is designed for analysis of > 1 omics experiments")
  }
  else {
    for (i in seq_along(c(1:n_omics))) {
      if (is.null(colnames(experiments(MAE_object)[[i]]))) {
        stop("All experiments must have colnames (sample names)")
      }
    }
  }

  if(harmonize){
    MAE_object <-
      MultiAssayExperiment::mergeReplicates(intersectColumns(MAE_object))
  }
  
  # Extract data matrices
  extracted_data <- MultiAssayExperiment::assays(MAE_object)@listData

  # Match experiment sample names to primary names if colData available
  if (length(colData(MAE_object)) > 0) {
    primary_names <- rownames(colData(MAE_object))

    for (i in seq_along(extracted_data)) {
      colnames(extracted_data[[i]]) <- primary_names
    }
  }
  rm(MAE_object)

  # Transpose to sample x feature format for input into MCIA
  extracted_data <- lapply(extracted_data, t)

  return(extracted_data)
}
