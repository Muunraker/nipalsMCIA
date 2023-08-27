#' Extract a list of harmonized data matrices from an MAE object
#'
#' @description Extract a list of harmonized data matrices for 
#' input into nipals_multiblock() from an MAE object
#' 
#' @param MAE_object an MAE object containing experiment data for extraction
#' colData field optional
#' experiments should either be SummarizedExperiment or SingleCellExperiment 
#' objects
#' @param subset_data \itemize{
#' \item `all` use all experiments in MAE object
#' \item `c(omic1,omic2,...)` list of omics from names(MAE_object)
#' }
#' @return List of harmonized data matrices for input into nipals_multiblock()
#' @examples
#' data(NCI60)
#' data_blocks_mae <- simple_mae(data_blocks,row_format="sample",
#'                               colData=metadata_NCI60)
#' NCI60_input = extract_from_mae(data_blocks_mae,subset='all')
#' @importFrom MultiAssayExperiment colData assays
#' @importFrom MultiAssayExperiment mergeReplicates intersectColumns
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @export


extract_from_mae <- function(MAE_object, subset_data="all"){

if (subset_data!="all"){
  if (intersect(subset_data,names(MAE_object))<2){
    stop("Please provide appropriate subset_data list")
  }
  else{
    MAE_object = MAE_object[,,subset]
  }
    
}

MAE_object_harmonize = MultiAssayExperiment::mergeReplicates(intersectColumns(MAE_object))

# Extract data matrices 
extracted_data <- MultiAssayExperiment::assays(MAE_object_harmonize)@listData 

# Match experiment sample names to primary names if colData available
if (length(colData(MAE_object_harmonize))>0){
  primary_names = rownames(colData(MAE_object_harmonize))
  for (i in 1:length(extracted_data)){colnames(extracted_data[[i]])=primary_names} 
}

# Transpose to sample x feature format for input into MCIA
extracted_data <- lapply(extracted_data, t) 

return(extracted_data)
}