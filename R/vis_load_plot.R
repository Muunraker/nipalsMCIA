#' Visualize all loadings on two factor axes
#'
#' @description Visualize all loadings recovered from
#' nipalsMCIA() output loadings matrix ranked using across two factor axes
#'
#' @param mcia_results object returned from nipals_multiblock() function
#' @param axes list of two numbers associated with two factors to visualize
#' @param color_pal a list of colors or function which returns a list of colors
#' @param color_pal_params a list of parameters for the color function
#' @return Plot of MCIA feature loadings for chosen axes
#' @examples
#' data(NCI60)
#' data_blocks_mae <- simple_mae(data_blocks,row_format="sample",
#'                               colData=metadata_NCI60)
#' mcia_results <- nipals_multiblock(data_blocks_mae, num_PCs = 10,
#'                                   plots = "none", tol = 1e-12)
#' vis_load_plot(mcia_results, axes = c(1, 4))
#' @importFrom ggplot2 ggplot aes_string geom_point theme_bw scale_color_manual
#' @importFrom ggplot2 xlab ylab labs
#' @export
vis_load_plot <- function(mcia_results, axes = c(1, 2),
                          color_pal = scales::viridis_pal,
                          color_pal_params = list()) {

    # get the colors
    plot_colors <- get_colors(mcia_results,
                              color_pal = color_pal,
                              color_pal_params = color_pal_params)


    # extracting global loadings
    gl <- mcia_results@global_loadings

    omic_dims <- vapply(mcia_results@block_loadings, dim, numeric(2))[1, ]
    omic_type <- c()
    omics_labels <- names(mcia_results@block_loadings)
    for (i in seq_along(omics_labels)) {
        omic_label <- omics_labels[i]
        length_omic <- omic_dims[i]
        omic_type <- c(omic_type, rep(omic_label, each = length_omic))
    }

    gl_f <- data.frame(gl[, axes])
    gl_f$omic <- omic_type
    gl_f$omic <- as.factor(gl_f$omic)
    colnames(gl_f) <- c(paste0("Axis_", axes[1]),
                        paste0("Axis_", axes[2]), "omic")

    # plot data
    omic <- gl_f$omic
    ggplot(data = gl_f,
           aes(x = .data[[colnames(gl_f)[1]]], y = .data[[colnames(gl_f)[2]]],
               color = omic)) +
      geom_point(alpha = 0.3) +
      labs(x = paste0("Axis ", axes[1]), y = paste0("Axis ", axes[2])) +
         # color = "Omic"
      scale_color_manual(values = plot_colors) +
      theme_bw()
}
