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
#' data(NCI60)
#' mcia_results <- nipals_multiblock(data_blocks, metadata = metadata_NCI60,
#' num_PCs = 10, plots = "none", tol = 1e-12)
#' colors_omics <- get_colors(mcia_results)
#' vis_load_plot(mcia_results, axes = c(1, 4), colors_omics = colors_omics)
#' @importFrom ggplot2 ggplot aes_string geom_point theme_bw scale_color_manual
#' @importFrom ggplot2 xlab ylab labs
#' @export

vis_load_plot <- function(mcia_out, axes = c(1, 2), colors_omics) {
  # extracting global loadings
  gl <- mcia_out$global_loadings

  # parsing data
  #omic_name <- gsub("^.*_", "", rownames(gl))
  
  omic_dims <- vapply(mcia_out$block_loadings,dim,numeric(2))[1,]
  omic_type<-c()
  omics_labels<-names(mcia_out$block_loadings)
  for (i in seq_along(omics_labels)){
    omic_label = omics_labels[i]
    length_omic = omic_dims[i]
    omic_type<-c(omic_type,rep(omic_label,each=length_omic))
  }
  
  gl_f <- data.frame(gl[, axes])
  gl_f$omic <- omic_type
  gl_f$omic <- as.factor(gl_f$omic)
  colnames(gl_f) <- c(paste0("Axis_", axes[1]),
                      paste0("Axis_", axes[2]),
                      "omic")

  # plot data
  p <- ggplot(data = gl_f,
              aes_string(x = colnames(gl_f)[1], y = colnames(gl_f)[2],
                         color = "omic")) +
         geom_point(alpha = 0.3) +
         labs(x = paste0("Axis ", axes[1]), y = paste0("Axis ", axes[2])) +
         scale_color_manual(values = colors_omics) +
         theme_bw()

  return(p)
}
