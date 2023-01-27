## ----setup-merger, include=FALSE----------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
set.seed(42) # NIPALS starts with a random vector

## ---- echo=FALSE, include=FALSE-----------------------------------------------
library(nipalsMCIA)
library(ggplot2)
library(ComplexHeatmap)

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
set.seed(42) # NIPALS starts with a random vector

## ---- echo=FALSE, include=FALSE-----------------------------------------------
library(nipalsMCIA)
library(ggplot2)
library(ComplexHeatmap)

## -----------------------------------------------------------------------------
# load the dataset, uses the name data_blocks
data(NCI60)
data_blocks$mrna[1:5,1:3]

data_blocks$miRNA[1:5,1:3]

data_blocks$prot[1:5,1:3]


## ---- warning=FALSE, message=FALSE--------------------------------------------
mcia_results <- nipals_multiblock(data_blocks, preprocMethod = "colprofile",
                                  num_PCs = 10, tol = 1e-12, plots = "none")

## -----------------------------------------------------------------------------
names(mcia_results)

## ----setup-V1-P1-I1, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
set.seed(42) # NIPALS starts with a random vector

## ---- echo=FALSE, include=FALSE-----------------------------------------------
library(nipalsMCIA)
library(ggplot2)
library(ComplexHeatmap)

## ---- echo=FALSE, include=FALSE-----------------------------------------------
if (!exists("data_blocks")){
    data(NCI60)
}

if (!exists("mcia_results")){
    mcia_results <- nipals_multiblock(data_blocks, preprocMethod='colprofile',
                                  num_PCs = 10, tol=1e-12, plots = 'none')
}

## ---- warning=FALSE, message=FALSE, fig.height=4, fig.width=8, fig.align='center'----
mcia_results <- nipals_multiblock(data_blocks, preprocMethod = "colprofile",
                                  num_PCs = 10, tol = 1e-12)

## ----Visualize Clusters, fig.height=4, fig.width=4, fig.align='center'--------
MCIA_plots(mcia_results, "projection_global", orders = c(1,2), legend_loc = "bottomleft")

## -----------------------------------------------------------------------------
metadata_NCI60

## ---- warning=FALSE, message=FALSE--------------------------------------------
# Conversion of metadata to data frame 
metadata_NCI60 <- data.frame(cancerType = unlist(metadata_NCI60))
row.names(metadata_NCI60) <- rownames(data_blocks[[1]])

# Adding metadata as part of the nipals_multiblock() function
mcia_results <- nipals_multiblock(data_blocks, preprocMethod = "colprofile", 
                                  metadata = metadata_NCI60, plots = "none",
                                  num_PCs = 10, tol = 1e-12)

# Alternative method for adding metadata
mcia_results$metadata <- metadata_NCI60

## ----Visualize Clusters Coloring, fig.height=4, fig.width=4, fig.align='center'----
MCIA_plots(mcia_results, "projection_global", orders = c(1,2), coloring = "cancerType", legend_loc = "bottomleft")

## ----setup-V1-P1-I2, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
set.seed(42) # NIPALS starts with a random vector

## ---- echo=FALSE, include=FALSE-----------------------------------------------
library(nipalsMCIA)
library(ggplot2)
library(ComplexHeatmap)

## ---- echo=FALSE, include=FALSE-----------------------------------------------
if (!exists("data_blocks")){
    data(NCI60)
}

if (!exists("mcia_results")){
    mcia_results <- nipals_multiblock(data_blocks, preprocMethod='colprofile',
                                  num_PCs = 10, tol=1e-12, plots = 'none')
}

