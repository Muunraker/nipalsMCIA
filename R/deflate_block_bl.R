#' Deflation via block loadings
#'
#' @description Removes data from a data frame in the direction of a given
#' block loadings vector.
#'
#' @details Subtracts the component of each row in the direction of a given
#' block loadings vector to yield a `deflated' data matrix.
#'
#' @param df a data frame in "sample" x "variable" format
#' @param bl a block loadings vector in variable space
#' @return the deflated data frame
#' @examples
#' df <- matrix(rbinom(15, 1, prob = 0.3), ncol = 3)
#' block_loading <- rbinom(3, 1, prob = 0.3)
#' deflated_data <- deflate_block_bl(df, block_loading)
#'
#' @export
deflate_block_bl <- function(df, bl) {
  df <- df - tcrossprod(as.matrix(df) %*% bl, bl)
  return(df)
}
