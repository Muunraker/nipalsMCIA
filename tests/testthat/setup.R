library(stringr)

# load in "data_blocks" and "metadata_NCI60"
data(NCI60)

# extract just the proteomics dataset
prot <- data_blocks$prot
# prot <- as.matrix(prot) # standardized proteomics dataset

# run MCIA with different numbers of scores
mcia_results <- nipals_multiblock(data_blocks, metadata = metadata_NCI60,
                                  col_preproc_method = "colprofile",
                                  num_PCs = 2, tol = 1e-12, plots = "none")

mcia_results3 <- nipals_multiblock(data_blocks, metadata = metadata_NCI60,
                                   col_preproc_method = "colprofile",
                                   num_PCs = 3, tol = 1e-12, plots = "none")

# run MCIA without metadata, used for testing the graphing section
mcia_results_no_meta <- nipals_multiblock(data_blocks,
                                          col_preproc_method = "colprofile",
                                          num_PCs = 2, tol = 1e-12,
                                          plots = "none")

# creating a shared set of meta colors
meta_colors <- get_metadata_colors(mcia_results, color_col = 1,
                                   color_pal_params = list(option = "E"))

# creating a shared set of omics colors
colors_omics <- get_colors(mcia_results)
