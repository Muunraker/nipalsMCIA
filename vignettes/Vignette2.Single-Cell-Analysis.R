## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----installation-github, eval = FALSE----------------------------------------
#  # devel version
#  
#  # install.packages("devtools")
#  devtools::install_github("Muunraker/nipalsMCIA", ref = "devel",
#                           force = TRUE, build_vignettes = TRUE) # devel version

## ----installation-bioconductor, eval = FALSE----------------------------------
#  # release version
#  if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#  
#  BiocManager::install("nipalsMCIA")

## ----load-packages, message = FALSE-------------------------------------------
# note that the TENxPBMCData package is not included in this list as you may
# decide to pull data from another source or use our provided objects

library(dplyr)
library(ggplot2)
library(ggpubr)
library(nipalsMCIA)
library(piggyback)
library(Seurat)

# NIPALS starts with a random vector
set.seed(42)

## ----set-paths----------------------------------------------------------------
# if you would like to save any of the data loaded/created locally
path_data <- file.path("..", "data")

# recommended location for external data
path_inst <- tempdir()

## ----pipelines, echo = FALSE--------------------------------------------------
# specify `tag = ` to use a different release other than latest

# pipeline image
piggyback::pb_download(file = "Vignette-2-Pipeline.png",
                       dest = path_inst, repo = "Muunraker/nipalsMCIA")

knitr::include_graphics(path = file.path(path_inst, "Vignette-2-Pipeline.png"))

## ----data-all-list------------------------------------------------------------
# list all of the currently available files in the latest release
piggyback::pb_list(repo = "Muunraker/nipalsMCIA", tag = "latest")

## ----data-all-download--------------------------------------------------------
# specify `tag = ` to use a different release other than latest

# files needed for running MCIA
piggyback::pb_download(file = "metadata_sc.csv",
                       dest = path_inst, repo = "Muunraker/nipalsMCIA")
# piggyback::pb_download(file = "data_blocks_sc.Rda",
#                        dest = path_inst, repo = "Muunraker/nipalsMCIA")

# MCIA results
piggyback::pb_download(file = "mcia_results_sc.Rds",
                       dest = path_inst, repo = "Muunraker/nipalsMCIA")

# marker genes for cell type annotation with Seurat
piggyback::pb_download(file = "marker_genes.csv",
                       dest = path_inst, repo = "Muunraker/nipalsMCIA")

# the Seurat data's metric summary file from 10x Genomics
piggyback::pb_download(file = "5k_pbmc_protein_v3_metrics_summary.csv",
                       dest = path_inst, repo = "Muunraker/nipalsMCIA")

# Seurat objects in different stages of processing
piggyback::pb_download(file = "tenx_pbmc5k_CITEseq_raw.rds",
                       dest = path_inst, repo = "Muunraker/nipalsMCIA")
piggyback::pb_download(file = "tenx_pbmc5k_CITEseq_annotated.rds",
                       dest = path_inst, repo = "Muunraker/nipalsMCIA")

## ----data-bioconductor-load, eval = FALSE-------------------------------------
#  # read in the data as a SingleCellExperiment object
#  tenx_pbmc3k <- TENxPBMCData::TENxPBMCData(dataset = "pbmc5k-CITEseq")
#  
#  # examine the data
#  tenx_pbmc3k
#  ## class: SingleCellExperiment
#  ## dim: 33538 5247
#  ## metadata(0):
#  ## assays(1): counts
#  ## rownames(33538): ENSG00000243485 ENSG00000237613 ... ENSG00000277475 ENSG00000268674
#  ## rowData names(4): ENSEMBL_ID Symbol_TENx Type Symbol
#  ## colnames: NULL
#  ## colData names(11): Sample Barcode ... Individual Date_published
#  ## reducedDimNames(0):
#  ## mainExpName: Gene Expression
#  ## altExpNames(1): Antibody Capture
#  
#  counts(tenx_pbmc3k)
#  ## <33538 x 5247> sparse matrix of class DelayedMatrix and type "integer":
#  ##                    [, 1]    [, 2]    [, 3]    [, 4] ... [, 5244] [, 5245] [, 5246] [, 5247]
#  ## ENSG00000243485       0       0       0       0   .       0       0       0       0
#  ## ENSG00000237613       0       0       0       0   .       0       0       0       0
#  ## ENSG00000186092       0       0       0       0   .       0       0       0       0
#  ## ENSG00000238009       0       0       0       0   .       0       0       0       0
#  ## ENSG00000239945       0       0       0       0   .       0       0       0       0
#  ##             ...       .       .       .       .   .       .       .       .       .
#  ## ENSG00000277856       0       0       0       0   .       0       0       0       0
#  ## ENSG00000275063       0       0       0       0   .       0       0       0       0
#  ## ENSG00000271254       0       0       0       0   .       0       0       0       0
#  ## ENSG00000277475       0       0       0       0   .       0       0       0       0
#  ## ENSG00000268674       0       0       0       0   .       0       0       0       0
#  
#  counts(altExp(tenx_pbmc3k))
#  ## <32 x 5247> sparse matrix of class DelayedMatrix and type "integer":
#  ##           [, 1]    [, 2]    [, 3]    [, 4] ... [, 5244] [, 5245] [, 5246] [, 5247]
#  ##    CD3      25     959     942     802   .     402     401       6    1773
#  ##    CD4     164     720    1647    1666   .    1417       1      46    1903
#  ##   CD8a      16       8      21       5   .       8     222       3       9
#  ##  CD11b    3011      12      11      11   .      15       7    1027       9
#  ##   CD14     696      12      13       9   .       9      17     382       8
#  ##    ...       .       .       .       .   .       .       .       .       .
#  ## HLA-DR     573      15      11      19   .       6      40     184      32
#  ##  TIGIT      10       3       3       3   .       2      15       1      12
#  ##   IgG1       4       4       2       4   .       1       0       2       4
#  ##  IgG2a       1       3       0       6   .       4       0       4       2
#  ##  IgG2b       6       2       4       8   .       0       0       2       5
#  
#  # examine the metadata:
#  head(colData(tenx_pbmc3k), n = 3)
#  ## DataFrame with 6 rows and 11 columns
#  ##           Sample            Barcode         Sequence   Library Cell_ranger_version Tissue_status Barcode_type
#  ##      <character>        <character>      <character> <integer>         <character>   <character>  <character>
#  ## 1 pbmc5k-CITEseq AAACCCAAGAGACAAG-1 AAACCCAAGAGACAAG         1              v3.0.2            NA     Chromium
#  ## 2 pbmc5k-CITEseq AAACCCAAGGCCTAGA-1 AAACCCAAGGCCTAGA         1              v3.0.2            NA     Chromium
#  ## 3 pbmc5k-CITEseq AAACCCAGTCGTGCCA-1 AAACCCAGTCGTGCCA         1              v3.0.2            NA     Chromium
#  ##     Chemistry Sequence_platform   Individual Date_published
#  ##   <character>       <character>  <character>    <character>
#  ## 1 Chromium_v3           NovaSeq HealthyDonor     2019-05-29
#  ## 2 Chromium_v3           NovaSeq HealthyDonor     2019-05-29
#  ## 3 Chromium_v3           NovaSeq HealthyDonor     2019-05-29
#  
#  head(rowData(tenx_pbmc3k), n = 3)
#  ## DataFrame with 6 rows and 4 columns
#  ##                      ENSEMBL_ID Symbol_TENx            Type       Symbol
#  ##                     <character> <character>     <character>  <character>
#  ## ENSG00000243485 ENSG00000243485 MIR1302-2HG Gene Expression           NA
#  ## ENSG00000237613 ENSG00000237613     FAM138A Gene Expression      FAM138A
#  ## ENSG00000186092 ENSG00000186092       OR4F5 Gene Expression        OR4F5
#  
#  metadata(tenx_pbmc3k)
#  ## list()
#  
#  # change the gene names from Ensembl IDs to the 10x genes
#  rownames(tenx_pbmc3k) <- rowData(tenx_pbmc3k)$Symbol_TENx

## ----data-bioconductor-formatting, eval = FALSE-------------------------------
#  # set up the list
#  data_blocks_sc_sce <- list()
#  data_blocks_sc_sce$mrna <- data.frame(as.matrix(counts(tenx_pbmc3k)))
#  data_blocks_sc_sce$adt <- data.frame(as.matrix(counts(altExp(tenx_pbmc3k))))
#  
#  summary(data_blocks_sc_sce)
#  ##      Length Class      Mode
#  ## mrna 5247   data.frame list
#  ## adt  5247   data.frame list
#  
#  # convert to a Seurat object (using `as.Seurat` won't work here)
#  obj_sce <- CreateSeuratObject(counts = data_blocks_sc_sce$mrna, # assay = "RNA"
#                                project = "pbmc5k_CITEseq")
#  obj_sce[["ADT"]] <- CreateAssayObject(counts = data_blocks_sc_sce$adt)
#  
#  # name the cells with their barcodes
#  obj_sce <- RenameCells(object = obj_sce,
#                         new.names = colData(tenx_pbmc3k)$Sequence)
#  
#  # add metadata from the SingleCellExperiment object
#  obj_sce <- AddMetaData(object = obj_sce,
#                         metadata = as.data.frame(colData(tenx_pbmc3k),
#                                                  row.names = Cells(obj_sce)))
#  
#  # this object will be slightly different than from the Seurat one down below
#  # e.g. 5297 rows vs. 4193 rows (since QC wasn't done) and different metadata
#  
#  head(obj_sce[[]], n = 3)
#  ##                     orig.ident nCount_RNA nFeature_RNA nCount_ADT nFeature_ADT         Sample            Barcode
#  ## AAACCCAAGAGACAAG SeuratProject       7375         2363       5178           31 pbmc5k-CITEseq AAACCCAAGAGACAAG-1
#  ## AAACCCAAGGCCTAGA SeuratProject       3772         1259       2893           29 pbmc5k-CITEseq AAACCCAAGGCCTAGA-1
#  ## AAACCCAGTCGTGCCA SeuratProject       4902         1578       3635           29 pbmc5k-CITEseq AAACCCAGTCGTGCCA-1
#  ##                          Sequence Library Cell_ranger_version Tissue_status Barcode_type   Chemistry Sequence_platform
#  ## AAACCCAAGAGACAAG AAACCCAAGAGACAAG       1              v3.0.2          <NA>     Chromium Chromium_v3           NovaSeq
#  ## AAACCCAAGGCCTAGA AAACCCAAGGCCTAGA       1              v3.0.2          <NA>     Chromium Chromium_v3           NovaSeq
#  ## AAACCCAGTCGTGCCA AAACCCAGTCGTGCCA       1              v3.0.2          <NA>     Chromium Chromium_v3           NovaSeq
#  ##                    Individual Date_published
#  ## AAACCCAAGAGACAAG HealthyDonor     2019-05-29
#  ## AAACCCAAGGCCTAGA HealthyDonor     2019-05-29
#  ## AAACCCAGTCGTGCCA HealthyDonor     2019-05-29
#  
#  # save the data locally if desired
#  save(data_blocks_sc_sce, obj_sce,
#       file = file.path(path_data, "data_sc_sce.Rda"))

## ----data-10x-object, eval = FALSE--------------------------------------------
#  # load the data (change the file path as needed)
#  data <- Seurat::Read10X(data.dir = file.path(path_data, "tenx_pbmc5k_CITEseq",
#                                               "filtered_feature_bc_matrix"),
#                          strip.suffix = TRUE) # remove the "-1"s from barcodes
#  ## 10X data contains more than one type and is being returned as a list
#  ## containing matrices of each type.
#  
#  # set minimum cells and/or features here if you'd like
#  obj <- Seurat::CreateSeuratObject(counts = data$`Gene Expression`,
#                                    project = "pbmc5k_CITEseq")
#  obj[["ADT"]] <- Seurat::CreateAssayObject(counts = data$`Antibody Capture`)
#  ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#  
#  # check the assays
#  Seurat::Assays(object = obj)
#  ## "RNA" "ADT"
#  
#  # list out the CITE-Seq surface protein markers
#  rownames(obj[["ADT"]])
#  ## [1] "CD3-TotalSeqB"            "CD4-TotalSeqB"           "CD8a-TotalSeqB"
#  ## [4] "CD11b-TotalSeqB"          "CD14-TotalSeqB"          "CD15-TotalSeqB"
#  ## [7] "CD16-TotalSeqB"           "CD19-TotalSeqB"          "CD20-TotalSeqB"
#  ## [10] "CD25-TotalSeqB"          "CD27-TotalSeqB"          "CD28-TotalSeqB"
#  ## [13] "CD34-TotalSeqB"          "CD45RA-TotalSeqB"        "CD45RO-TotalSeqB"
#  ## [16] "CD56-TotalSeqB"          "CD62L-TotalSeqB"         "CD69-TotalSeqB"
#  ## [19] "CD80-TotalSeqB"          "CD86-TotalSeqB"          "CD127-TotalSeqB"
#  ## [22] "CD137-TotalSeqB"         "CD197-TotalSeqB"         "CD274-TotalSeqB"
#  ## [25] "CD278-TotalSeqB"         "CD335-TotalSeqB"         "PD-1-TotalSeqB"
#  ## [28] "HLA-DR-TotalSeqB"        "TIGIT-TotalSeqB"         "IgG1-control-TotalSeqB"
#  ## [31] "IgG2a-control-TotalSeqB" "IgG2b-control-TotalSeqB"
#  
#  # save the data locally if desired
#  saveRDS(obj, file.path(path_data, "tenx_pbmc5k_CITEseq_raw.rds"))

## ----mcia-metadata------------------------------------------------------------
# read in the annotated cells
metadata_sc <- read.csv(file = file.path(path_inst, "metadata_sc.csv"),
                        header = TRUE, row.names = 1)

# examples
metadata_sc %>% slice_sample(n = 5)

## ----mcia-decomp-load-data, eval = FALSE--------------------------------------
#  # load the object setup for running MCIA [10x Genomics & Seurat]
#  load(file = file.path(path_inst, "data_blocks_sc.Rda"))

## ----mcia-decomp-run, eval = FALSE--------------------------------------------
#  # "largest_sv" results in a more balanced contribution
#  # from the blocks than the default "unit_var"
#  set.seed(42)
#  
#  # convert data_blocks_sc to an MAE object using the SingleCellExperiment class
#  data_blocks_sc_mae <-
#    MultiAssayExperiment::MultiAssayExperiment(lapply(data_blocks_sc, function(x)
#      SingleCellExperiment::SingleCellExperiment(t(as.matrix(x)))),
#      colData = metadata_sc)
#  mcia_results_sc <- nipals_multiblock(data_blocks = data_blocks_sc_mae,
#                                       col_preproc_method = "colprofile",
#                                       block_preproc_method = "largest_sv",
#                                       num_PCs = 10, tol = 1e-9,
#                                       deflationMethod = "global",
#                                       plots = "none")
#  ## Performing column-level pre-processing...
#  ## Column pre-processing completed.
#  ## Performing block-level preprocessing...
#  ## Block pre-processing completed.
#  ## Computing order 1 scores
#  ## Computing order 2 scores
#  ## Computing order 3 scores
#  ## Computing order 4 scores
#  ## Computing order 5 scores
#  ## Computing order 6 scores
#  ## Computing order 7 scores
#  ## Computing order 8 scores
#  ## Computing order 9 scores
#  ## Computing order 10 scores
#  
#  # saveRDS(mcia_results_sc, file = file.path(path_data, "mcia_results_sc.Rds"))

## ----mcia-decomp-load---------------------------------------------------------
# load the results of the previous block (if already run and saved)
mcia_results_sc <- readRDS(file = file.path(path_inst, "mcia_results_sc.Rds"))
mcia_results_sc

## ----mcia-plots-colors--------------------------------------------------------
# for the projection plot
# technically you could just do color_pal_params = list(option = "D"), but saving
# the colors is useful for other plots like in the Seurat section
meta_colors_sc <- get_metadata_colors(mcia_results = mcia_results_sc,
                                      color_col = "CellType",
                                      color_pal = scales::viridis_pal,
                                      color_pal_params = list(option = "D"))

# for other plots
colors_omics_sc <- get_colors(mcia_results = mcia_results_sc)

## ----mcia-plots-eigenvalue-scree, fig.dim = c(5, 4)---------------------------
global_scores_eigenvalues_plot(mcia_results = mcia_results_sc)

## ----mcia-plots-projection, fig.dim = c(5, 5)---------------------------------
projection_plot(mcia_results = mcia_results_sc,
                projection = "global", orders = c(1, 2),
                color_col = "CellType", color_pal = meta_colors_sc,
                legend_loc = "bottomright")

## ----mcia-plots-heatmap-global, fig.dim = c(7, 5)-----------------------------
suppressMessages(global_scores_heatmap(mcia_results = mcia_results_sc,
                                       color_col = "CellType",
                                       color_pal = meta_colors_sc))

## ----mcia-plots-heatmap-block, fig.dim = c(4, 2.5)----------------------------
block_weights_heatmap(mcia_results = mcia_results_sc)

