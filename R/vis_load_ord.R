#' Visualize ranked loadings
#'
#' @description Visualize a scree plot of loadings recovered from 
#' nipalsMCIA() output loadings matrix ranked using the 
#' ord_loadings() functions
#'
#' @param gl_f_ord Ranked loading dataframe output from ord_loadings() function
#' @param colors_omics named list of colors associated with omics,
#' output of get_colors() function
#' @param n_feat number of features to visualize
#' @examples
#' vis_load_ord(all_pos_1_df, colors_omics = colors_omics)
#' @importFrom ggplot2 ggplot
#' @export

vis_load_ord <- function(gl_f_ord, colors_omics, n_feat = 15) {
  n_plot <- min(nrow(gl_f_ord), n_feat)
  omic_subset <- names(table(droplevels(gl_f_ord[0:n_plot,])$omic))
  color_vals <- colors_omics[omic_subset]
  
  p <- ggplot2::ggplot(data = gl_f_ord[1:n_plot,], 
                       aes(x = factor(omic_name, level = omic_name), 
                           y = loading, color = omic)) +
    geom_point() +
    theme_bw() +
    xlab('Feature') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
    scale_color_manual(values = color_vals)
  
  return(p)
}