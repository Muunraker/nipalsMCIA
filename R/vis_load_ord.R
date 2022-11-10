#' Visualize ranked loadings
#'
#' @description Visualize a scree plot of loadings recovered from
#' nipalsMCIA() output loadings matrix ranked using the
#' ord_loadings() functions
#'
#' @param gl_f_ord Ranked loading dataframe output from ord_loadings() function
#' @omic name of the given omic dataset
#' @param colors_omics named list of colors associated with omics,
#' output of get_colors() function
#' @param n_feat number of features to visualize
#' @return Plot in features for a factor by rank
#' @examples
#' vis_load_ord(all_pos_1_df, colors_omics = colors_omics)
#' @importFrom ggplot2 ggplot aes geom_point theme_bw xlab theme
#' @importFrom ggplot2 scale_color_manual
#' @export

vis_load_ord <- function(gl_f_ord, omic_name, colors_omics, n_feat = 15) {
  loading <- rlang::quo(`loading`)
  omic <- rlang::quo(`omic`)
  n_plot <- min(nrow(gl_f_ord), n_feat)
  omic_subset <- names(table(droplevels(gl_f_ord[seq(0, n_plot), ])$omic))
  color_vals <- colors_omics[omic_subset]
  p <- ggplot(data = gl_f_ord[seq_len(n_plot), ],
              ggplot2::aes(x = factor(omic_name, level = omic_name),
                  y = !!loading, color = !!omic)) +
       ggplot2::geom_point() +
       ggplot2::theme_bw() +
       ggplot2::xlab("Feature") +
       ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                 hjust = 1,
                                                 size = 6)) +
      ggplot2::scale_color_manual(values = color_vals)
  return(p)
}
