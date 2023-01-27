#' Perform biological annotation-based comparison
#'
#' @description Runs fgsea for the input gene vector
#' @param metagenes Vector of gene scores where the row names are HUGO symbols
#' @param path.database path to a GMT annotation file
#' @param factors vector of factors which should be analyzed
#' @param pval.thr p-value threshold (default to 0.05)
#' @param nproc number of processors to utilize
#' @return data frame with the most significant p-value number of significant
#' pathways
#' @return the selectivity scores across the given factors
#' @importFrom stats setNames
#' @export
gsea_report <- function(metagenes, path.database, factors = NULL,
                        pval.thr = 0.05, nproc = 4) {
  # Load annotation database
  pathways <- fgsea::gmtPathways(path.database)

  # Number of factors
  if (is.null(factors)) {
    factors <- seq_len(ncol(metagenes))
  }

  # Containers to report results
  report_min_pval <- numeric(0)
  report_number <- numeric(0)
  sig_paths <- numeric(0)

  # Calculate biological annotation enrichment, for each factor
  for (j in factors) {
    print(paste0("Running GSEA for Factor", j))

    # Extract a vector of scores for GSEA and set rownams
    # to HUGO symbols
    scores <- setNames(as.matrix(metagenes[, j]), rownames(metagenes))

    # Compute GSEA
    fgseaRes <- fgsea::fgseaMultilevel(pathways, scores,
                                       nPermSimple = 10000, minSize = 15,
                                       maxSize = 500, nproc = nproc)

    # Report if at least one pathway is significant and
    # the min-pvalue from all pathways
    if (sum(fgseaRes$padj < pval.thr) != 0) {
      # Report the minimum p-value
      min_pval <- min(fgseaRes$padj)
      report_min_pval <- rbind(report_min_pval, min_pval)

      # Keep names of significant pathways
      curr_sig_paths <- fgseaRes[fgseaRes$padj < pval.thr, "pathway"]
      sig_paths <- c(sig_paths, curr_sig_paths)

      # Report number of unique significant pathways
      report_number <- rbind(report_number, dplyr::n_distinct(curr_sig_paths))

    }
    else {
      # Report the minimum p-value, assigning NA
      report_min_pval <- rbind(report_min_pval, NA)

      # Report number of unique significant pathways, assigning 0
      report_number <- rbind(report_number, 0)
    }
  }

  # Report selectivity
  # Selectivity is calculated best across ALL factors
  # The formula is give within https://go.nature.com/3EjW5HJ
  # and section "Selectivity Score". SS = (Nc + Nf) / 2L
  # Nc is the total number of clinical annotations associated with at
  # least a factor, Nf the total number of factors associated with at
  # least a clinical annotation, and L the total number of
  # associations between clinical annotations and factors. S has a
  # maximum value of 1 when each factor is associated with one and
  # only one clinical/biological annotation, and a minimum of 0 in
  # the opposite case. An optimal method should thus maximize its
  # number of factors associated with clinical/biological annotations
  # without having a too low selectivity.
  Nc <- dplyr::n_distinct(sig_paths)
  Nf <- length(which(!is.na(min_pval)))
  L <- length(sig_paths)
  SS <- (Nc + Nf) / (2 * L)

  # Setting up the final reporting data structure with list containing
  # two values, the first is a dataframe where the rows are factors and
  # the columns is composed of the following fields:
  # - the minimum p-value of all pathways
  # - the total number of significant pathways
  # Lastly, the second value is the selectivity value which assesses whether
  # the set of factors are capturing very different gene sets/pathways from
  # one another
  out <- list()
  out[[1]] <- data.frame(min_pval = report_min_pval,
               total_pathways = report_number)
  row.names(out[[1]]) <- paste0("Factor", factors)

  out[[2]] <- SS
  names(out) <- c("per-factor-results", "selectivity")
  return(out)
}
