data(NCI60)
mcia_results <- nipals_multiblock(data_blocks, metadata = metadata_NCI60,
                                  preproc_method = "colprofile",
                                  num_PCs = 2, tol = 1e-12, plots = "none")
