library(nipalsMCIA)
source("correlate_factors_to_features.R")
data(NCI60)

######## Run basic MCIA example ########
results <- nipals_multiblock(data_blocks, preprocMethod='colprofile',num_PCs = 10,
                             plots='global', tol=1e-12)
feature_mat <- data_blocks[["mrna"]]
latent_mat <- results$global_scores

View(data_blocks)
View(results)


######## Test corr_heatmap_fvl_ggplot2() ########
corr_heatmap_fvl_ggplot2(latent_mat, feature_mat[,0:20])

######## Test corr_heatmap_fvl_ComplexHeatmap() ########
corr_heatmap_fvl_ComplexHeatmap(latent_mat, feature_mat[,0:20])

######## global_loadings_heatmap_ComplexHeatmap() ########
global_loadings_heatmap_ComplexHeatmap(results$global_loadings,
                                       data_blocks,
                                       "mrna",
                                       seq(1,10))
    
######## global_loadings_heatmap_ComplexHeatmap() ########
global_loadings_heatmap_ComplexHeatmap(results$global_loadings,
                                       data_blocks,
                                       "miRNA",
                                       seq(1,10))

######## global_loadings_heatmap_ComplexHeatmap() ########
global_loadings_heatmap_ComplexHeatmap(results$global_loadings,
                                       data_blocks,
                                       "prot",
                                       seq(1,10))


######## plot global scores ########
gs_scores = results$global_scores
colnames(gs_scores) = paste0('LF', seq(1, ncol(results$global_scores)))
p = Heatmap(gs_scores, 
    name = "GS Score", 
    column_title = "Latent Factors",
    row_title = "Samples",
    row_names_gp = gpar(fontsize = 7),
    show_column_names = T,
    show_row_names = T,
    row_names_side = "right"
)    
p
