#' Centered Column Profile Pre-processing
#'
#' @description Converts data blocks into centered column profiles where each
#' block has unit variance. Mimics the pre-processing in the Omicade4 package
#' (Meng et al. 2014)
#' @details Performs the following steps on a given data frame: \itemize{
#' \item Offsets data to make whole matrix non-negative
#' \item Divides each column by its sum
#' \item Subtracts (row sum/total sum) from each row
#' \item Multiplies each column by sqrt(column sum/total sum)
#' \item Divides the whole data frame by its total variance (the sqrt of the
#' sum of singular values)
#' }
#' @param df the data frame to apply pre-processing to, in "sample" x
#' "variable" format
#' @return the processed data frame
#' @examples
#' preprocessed_dataframe <- CCpreproc(dataframe)
#' @export
CCpreproc <- function(df) {
  temp_df <- as.matrix(df)

  # Making data non-negative
  minVal <- min(temp_df)
  if (minVal < 0) {
    offset <- floor(minVal)
    temp_df <- temp_df + abs(offset)
  }

  # Generating centered column profiles:
  totsum <- sum(temp_df)
  colsums <- colSums(temp_df)
  row_contribs <- rowSums(temp_df) / totsum

  ## Dividing by column sums
  nz_cols <- which(colsums != 0) # excluding zero columns to avoid NaNs
  temp_df[, nz_cols] <- t(t(temp_df[, nz_cols]) / colsums[nz_cols])

  ## Subtracting row contributions
  temp_df <- temp_df - row_contribs

  # Applying feature weighting ("multiplication by feature metrics")
  temp_df <- t(t(temp_df) * sqrt(colsums / totsum))

  # Weighting block to unit variance (block variance computed via Fro norm)
  block_var <- norm(temp_df, type = "F")^2 / (max(1, nrow(df) - 1))

  if (block_var != 0) {
    temp_df <- temp_df * (1 / sqrt(block_var))
  } else {
    warning("New data has zero variance.")
  }
  return(temp_df)
}
