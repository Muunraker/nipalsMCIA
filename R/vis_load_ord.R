#' Visualize ranked loadings
#'
#' @description Visualize a scree plot of loadings recovered from
#' nipalsMCIA() output loadings matrix ranked using the
#' ord_loadings() functions
#'
#' @param gl_f_ord Ranked loading dataframe output from ord_loadings() function
#' @param omic_name name of the given omic dataset
#' @param colors_omics named list of colors associated with omics,
#' output of get_colors() function
#' @param n_feat number of features to visualize
#' @return Plot in features for a factor by rank
#' @examples
#' data(NCI60)
#' mcia_results <- nipals_multiblock(data_blocks, metadata = metadata_NCI60,
#' num_PCs = 10, plots = "none", tol = 1e-12)
#' all_pos_1 <- ord_loadings(mcia_out = mcia_results, omic = "all",
#'   absolute = FALSE, descending = TRUE, factor = 1)
#' colors_omics <- get_colors(mcia_results)
#' vis_load_ord(all_pos_1, colors_omics = colors_omics)
#' @importFrom rlang quo
#' @importFrom ggplot2 ggplot aes geom_point theme_bw xlab theme
#' @importFrom ggplot2 scale_color_manual element_text
#' @export

vis_load_ord <- function(gl_f_ord, omic_name, colors_omics, n_feat = 15) {
    
  # setting parameters with tidy evaluation
  loading <- quo(`loading`)
  omic <- quo(`omic`)
  
  # making a plot for the given omics
  n_plot <- min(nrow(gl_f_ord), n_feat)
  omic_subset <- names(table(droplevels(gl_f_ord[seq(0, n_plot), ])$omic))
  color_vals <- colors_omics[omic_subset]
  p <- ggplot(data = gl_f_ord[seq_len(n_plot), ],
              aes(x = factor(omic_name, level = omic_name),
                  y = !!loading, color = !!omic)) +
       geom_point() +
       theme_bw() +
       xlab("Feature") +
       theme(axis.text.x = element_text(angle = 45,
                                        hjust = 1,
                                        size = 6)) +
      scale_color_manual(values = color_vals)
  return(p)
}
