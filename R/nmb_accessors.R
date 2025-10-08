#' Accessor function for global scores
#'
#' @description Retrieves the global scores as a matrix from a `NipalsResult`
#' object, typically output from `nipals_multiblock()`.
#' @param nmb_object A `NipalsResult` object.
#' @return a matrix containing global scores.
#' @examples
#' data("NCI60")
#' data_blocks_mae <- simple_mae(data_blocks, row_format = "sample",
#'                               colData = metadata_NCI60)
#' mcia_out <- nipals_multiblock(data_blocks_mae, num_PCs = 10)
#' global_scores <- nmb_get_gs(mcia_out)
#' @export

nmb_get_gs <- function(nmb_object) {
    return(nmb_object@global_scores)
}


#' Accessor function for global loadings
#'
#' @description Retrieves the global loadings as a matrix from a `NipalsResult`
#' object, typically output from `nipals_multiblock()`.
#' @param nmb_object A `NipalsResult` object.
#' @return a matrix containing global loadings.
#' @examples
#' data("NCI60")
#' data_blocks_mae <- simple_mae(data_blocks, row_format = "sample",
#'                               colData = metadata_NCI60)
#' mcia_out <- nipals_multiblock(data_blocks_mae, num_PCs = 10)
#' global_loadings <- nmb_get_gl(mcia_out)
#' @export

nmb_get_gl <- function(nmb_object) {
    return(nmb_object@global_loadings)
}


#' Accessor function for block scores
#'
#' @description Retrieves the block scores as a list of matrices
#' from a `NipalsResult` object, typically output from `nipals_multiblock()`.
#' @param nmb_object A `NipalsResult` object.
#' @return a list of matrices containing block scores.
#' @examples
#' data("NCI60")
#' data_blocks_mae <- simple_mae(data_blocks, row_format = "sample",
#'                               colData = metadata_NCI60)
#' mcia_out <- nipals_multiblock(data_blocks_mae, num_PCs = 10)
#' block_scores <- nmb_get_bs(mcia_out)
#' @export

nmb_get_bs <- function(nmb_object) {
    return(nmb_object@block_scores)
}


#' Accessor function for block loadings
#'
#' @description Retrieves the block loadings as a list of matrices
#' from a `NipalsResult` object, typically output from `nipals_multiblock()`.
#' @param nmb_object A `NipalsResult` object.
#' @return a list of matrices containing block loadings.
#' @examples
#' data("NCI60")
#' data_blocks_mae <- simple_mae(data_blocks, row_format = "sample",
#'                               colData = metadata_NCI60)
#' mcia_out <- nipals_multiblock(data_blocks_mae, num_PCs = 10)
#' block_loadings<- nmb_get_bl(mcia_out)
#' @export

nmb_get_bl <- function(nmb_object) {
    return(nmb_object@block_loadings)
}


#' Accessor function for block score weights
#'
#' @description Retrieves the block score weights
#' from a `NipalsResult` object, typically output from `nipals_multiblock()`.
#' @param nmb_object A `NipalsResult` object.
#' @return a matrix containing the block score weights.
#' @examples
#' data("NCI60")
#' data_blocks_mae <- simple_mae(data_blocks, row_format = "sample",
#'                               colData = metadata_NCI60)
#' mcia_out <- nipals_multiblock(data_blocks_mae, num_PCs = 10)
#' block_score_weights <- nmb_get_bs_weights(mcia_out)
#' @export

nmb_get_bs_weights <- function(nmb_object) {
    return(nmb_object@block_score_weights)
}

#' Accessor function for eigenvalues
#'
#' @description Retrieves the eigenvalues
#' from a `NipalsResult` object, typically output from `nipals_multiblock()`.
#' @param nmb_object A `NipalsResult` object.
#' @return a matrix containing the eigenvalues for all global scores.
#' @examples
#' data("NCI60")
#' data_blocks_mae <- simple_mae(data_blocks, row_format = "sample",
#'                               colData = metadata_NCI60)
#' mcia_out <- nipals_multiblock(data_blocks_mae, num_PCs = 10)
#' nipals_eigvals <- nmb_get_eigs(mcia_out)
#' @export

nmb_get_eigs <- function(nmb_object) {
    return(nmb_object@eigvals)
}

#' Accessor function for metadata
#'
#' @description Retrieves the metadata
#' from a `NipalsResult` object, typically output from `nipals_multiblock()`.
#' @param nmb_object A `NipalsResult` object.
#' @return a dataframe containing metadata associated with the `NipalsResult` object.
#' @examples
#' data("NCI60")
#' data_blocks_mae <- simple_mae(data_blocks, row_format = "sample",
#'                               colData = metadata_NCI60)
#' mcia_out <- nipals_multiblock(data_blocks_mae, num_PCs = 10)
#' nipals_metadata <- nmb_get_metadata(mcia_out)
#' @export

nmb_get_metadata <- function(nmb_object) {
  return(nmb_object@metadata)
}
