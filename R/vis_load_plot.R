#' Visualize all loadings on two factor axes
#'
#' @description Visualize all loadings recovered from
#' nipalsMCIA() output loadings matrix ranked using across two factor axes
#'
#' @param mcia_out object returned from nipals_multiblock() function
#' @param axes list of two numbers associated with two factors to visualize
#' @param colors_omics named list of colors associated with omics,
#' output of get_colors() function
#' @return Plot of MCIA feature loadings for chosen axes
#' @examples
#' vis_load_plot(mcia_results, axes = c(1, 4), colors_omics)
#' @importFrom ggplot2 ggplot aes_string geom_point theme_bw scale_color_manual
#' @export

vis_load_plot <- function(mcia_out, axes = c(1, 2), colors_omics) {
  gl <- mcia_out$global_loadings
  omic_name <- gsub("^.*_", "", rownames(gl))
  gl_f <- data.frame(gl[, axes])
  gl_f$omic <- omic_name
  gl_f$omic <- as.factor(gl_f$omic)
  colnames(gl_f) <- c(paste0("Axis_", axes[1]),
                      paste0("Axis_", axes[2]),
                      "omic")
  p <- ggplot(data = gl_f,
              aes_string(x = colnames(gl_f)[1],
              y = colnames(gl_f)[2], color = "omic")) +
    geom_point(alpha = 0.3) +
    theme_bw() +
    scale_color_manual(values = colors_omics)
  return(p)
}
