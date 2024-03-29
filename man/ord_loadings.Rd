% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ord_loadings.R
\name{ord_loadings}
\alias{ord_loadings}
\title{Ranked global loadings dataframe}
\usage{
ord_loadings(
  mcia_results,
  omic = "all",
  factor = 1,
  absolute = FALSE,
  descending = TRUE
)
}
\arguments{
\item{mcia_results}{object returned from nipals_multiblock() function}

\item{omic}{choose an omic to rank, or choose 'all' for all,
((omic = "all", omic = "miRNA", etc.))}

\item{factor}{choose a factor (numeric value from 1 to number of factors
in mcia_results)}

\item{absolute}{whether to rank by absolute value}

\item{descending}{whether to rank in descending or ascending order}
}
\value{
ranked dataframe
}
\description{
Creates a dataframe with ranked loadings for a given factor
}
\examples{
data(NCI60)
data_blocks_mae <- simple_mae(data_blocks,row_format="sample",
                              colData=metadata_NCI60)
mcia_results <- nipals_multiblock(data_blocks_mae, num_PCs = 10,
                                  plots = "none", tol = 1e-12)
all_pos_1 <- ord_loadings(mcia_results = mcia_results, omic = "all",
                          absolute = FALSE, descending = TRUE, factor = 1)
}
