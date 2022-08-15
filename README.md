
---
title: "NIPALS-MCIA: Software to Compute Multi-Block Dimensionality Reduction via Nonlinear Iterative Partial Least Squares"
output:
  html_document:
    toc: true
---

## Installation Instructions 

### Required Dependencies

* The R component of NIPALS-MCIA requires R versions 3.5.0 or above. Some functions from the Pracma R package are automatically installed with NIPALS-MCIA.


### Installation

The R version of MIPALS-MCIA can be easily installed from Github using the `devtools' package. 

```{r}
install.packages("devtools")
library(devtools)
install_github("Muunraker/NIPALS-MCIA",force = TRUE)

```
The MATLAB version of NIPALS-MCIA is installed simply by downloading the files from Github 
and adding the 'Functions' folder to the MATLAB path.

```{matlab}
addpath('<path to MATLAB_MCIA>\Functions')

```

## Using NIPALS-MCIA

### Performing dimensionality-reduction

The package currently includes only one vignette - a cut version of a three-omic NCI-60
cancer cell line dataset. Only one function call is needed to perform MCIA: 
```{r}
data(NCI60) # import data as "data_blocks"
mcia_results <- nipals_multiblock(data_blocks, preprocMethod='colprofile',num_PCs = 8, tol=1e-9)

```
Note: this dataset is reproduced from the omicade4 package (Meng et. al., 2014). This package assumes
all input datasets are in sample by feature format. 

Meng, C., Kuster, B., Culhane, A. C., & Gholami, A. M. (2014). A multivariate approach to the integration of multi-omics datasets. BMC bioinformatics, 15(1), 1-13.

### Predicting scores on new data

NIPALS-MCIA includes the 'predict_gs' function, which allows for the computation of global
scores for new data given a previously-computed set of loadings.
```{r}
# Computing MCIA on a dataset:
mcia_results <- nipals_multiblock(data_blocks, preprocMethod='colprofile',num_PCs = 8, tol=1e-9)

# Prediction
bl <- mcia_results$block_loadings # existing block loadings matrix
bw <- mcia_results$block_score_weights  # existing block weights matrix
new_data <- list(omics_1, omics_2, omics_3) # new data where each omics is a matrix or dataframe.
gs_new <- predict_gs(bl, bw, new_data) # new global score from "new_data" 


```
Note that the new dataset must have the same number of omics and match all the features of the original dataset. 
The new dataset must also be formatted as a list of data frames or data matrices in sample by feature form.
