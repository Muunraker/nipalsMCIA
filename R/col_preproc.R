#' Centered Column Profile Pre-processing
#'
#' @description Converts data blocks into centered column profiles where each
#' block has unit variance. Mimics the pre-processing in the Omicade4 package
#' (Meng et al. 2014)
#' @details Performs preprocessing on a sample/variable (row/column) level
#' according to the parameter given. 
#' @param df the data frame to apply pre-processing to, in "sample" x
#' "variable" format
#' @param col_preproc_method denotes the type of column-centered preprocessing.
#' Options are: \itemize{
#' \item `colprofile` Performs the following steps on a given data frame: 
#' \itemize{
#' \item Offsets data to make whole matrix non-negative
#' \item Divides each column by its sum
#' \item Subtracts (row sum/total sum) from each row
#' \item Multiplies each column by sqrt(column sum/total sum)
#' }
#' \item `standardized` centers each column and divides by its standard 
#' deviation. 
#' \item `centered_only` ONLY centers data
#' }
#' @return the processed data frame
#' @examples
#' df <- matrix(rbinom(15, 1, prob = 0.3), ncol = 3) 
#' preprocessed_dataframe <- col_preproc(df, col_preproc_method = 'colprofile')
#' @export
col_preproc <- function(df, col_preproc_method) {
  temp_df <- as.matrix(df)
  
  if(tolower(col_preproc_method) == "colprofile"){
    
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
    
    # Dividing by column sums
    nz_cols <- which(colsums != 0) # excluding zero columns to avoid NaNs
    temp_df[, nz_cols] <- t(t(temp_df[, nz_cols]) / colsums[nz_cols])
    
    # Subtracting row contributions
    temp_df <- temp_df - row_contribs
    
    # Applying feature weighting by proportion of variance
    temp_df <- t(t(temp_df) * sqrt(colsums / totsum))
    
  }else if(tolower(col_preproc_method) == "standardized"){
    temp_df <- scale(temp_df)
      
  }else if(tolower(col_preproc_method) == "centered_only"){
    temp_df <- scale(temp_df ,center = TRUE, scale = FALSE)
    
  }else{
    stop("Column preprocessing method not recognized - pick from available options")
  }
  
  
  return(temp_df)
}
