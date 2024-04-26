#' Visualize ranked loadings
#'
#' @description Visualize a scree plot of loadings recovered from
#' nipalsMCIA() output loadings matrix ranked using the
#' ord_loadings() functions
#'
#' @param mcia_results object returned from nipals_multiblock() function
#' @param omic name of the given omic dataset
#' @param factor choose a factor (numeric value from 1 to number of factors
#'               in mcia_results)
#' @param n_feat number of features to visualize
#' @param absolute whether to rank by absolute value
#' @param descending whether to rank in descending or ascending order
#' @param color_pal a list of colors or function which returns a list of colors
#' @param color_pal_params a list of parameters for the color function
#' @return Plot in features for a factor by rank
#' @examples
#' data(NCI60)
#' data_blocks_mae <- simple_mae(data_blocks,row_format="sample",
#'                               colData=metadata_NCI60)
#' mcia_results <- nipals_multiblock(data_blocks_mae, num_PCs = 10,
#'                                   plots = "none", tol = 1e-12)
#' vis_load_ord(mcia_results, omic="mrna")
#' @importFrom rlang quo
#' @importFrom ggplot2 ggplot aes geom_point theme_bw xlab theme
#' @importFrom ggplot2 scale_color_manual element_text
#' @export
vis_load_ord <- function(mcia_results, omic, factor = 1, n_feat = 15,
                         absolute = TRUE, descending = TRUE,
                         color_pal = scales::viridis_pal,
                         color_pal_params = list()) {
    # get the colors
    plot_colors <- get_colors(mcia_results,
                              color_pal = color_pal,
                              color_pal_params = color_pal_params)

    # get the load ordering internally
    gl_f_ord <- ord_loadings(mcia_results = mcia_results, omic = omic,
                             factor = factor, absolute = absolute,
                             descending = descending)
    omic_name <- gl_f_ord$omic_name

    # setting parameters with tidy evaluation
    loading <- quo(`loading`)
    omic_quo <- quo(`omic`)

    # making a plot for the given omics
    n_plot <- min(nrow(gl_f_ord), n_feat)
    omic_subset <- names(table(droplevels(gl_f_ord[seq(0, n_plot), ])$omic))
    color_vals <- plot_colors[omic_subset]

    # plot the data
    ggplot(data = gl_f_ord[seq_len(n_plot), ],
           aes(x = factor(omic_name, levels = omic_name),
               y = !!loading, color = !!omic_quo)) +
        geom_point() +
        labs(x = "Feature", y = "Loading") + # color = "Omic"
        scale_color_manual(values = color_vals, limits = force) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
}
