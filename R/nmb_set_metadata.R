#' Function to set the metadata of a `NipalsResult` object.
#'
#' @param nmb_object The `NipalsResult` object to set metadata for.
#' @param nmb_metadata A data frame of metadata for each sample in the dataset.
#' Each column must have length equal to the number of samples (rows)
#' in the dataset.
#' @return A `NipalsResult` object with the provided metadata.
#' @examples
#' data("NCI60")
#' mcia_out <- nipals_multiblock(data_blocks, num_PCs = 10)
#' mcia_out <- nmb_set_metadata(mcia_out,metadata_NCI60)
#' @export

nmb_set_metadata <- function(nmb_object, nmb_metadata) {
    # check metadata has correct type
    if(!is(nmb_metadata,"data.frame")){
        errmsg <- sprintf("Metadata must be a data frame.")
        stop(errmsg)
    }

    # check metadata is nonempty
    if(any(dim(nmb_metadata)==0)){
        errmsg <- sprintf("Metadata must be nonempty.")
        stop(errmsg)
    }

    # check metadata has matching number of rows (Samples)
    if(nrow(nmb_metadata) != nrow(nmb_object@global_scores)){
        errmsg <- sprintf("Number of rows in metadata does not match dataset.")
        stop(errmsg)
    }

    # Check that metadata sample names match block sample names
    if (any(tolower(rownames(nmb_object@global_scores)) != tolower(rownames(nmb_metadata)))) {
        errmsg <- sprintf("Metadata sample names dont match dataset sample names.")
        warning(errmsg)
    }
    # set metadata
    nmb_object@metadata <- nmb_metadata

    return(nmb_object)
}
