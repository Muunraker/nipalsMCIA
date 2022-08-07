# NIPALS-MCIA: Software to Compute Multi-Block Dimensionality Reduction via Nonlinear Iterative Partial Least Squares 

[[_TOC_]]

# Installation Instructions 

### Required Dependencies

* The R component of NIPALS-MCIA requires R versions 3.5.0 or above. Some functions from the Pracma R package are automatically installed with NIPALS-MCIA.
* The MATLAB component of NIPALS-MCIA was built on version R2020b and later. 


### Installation

The R version of MIPALS-MCIA can be easily installed from Github using the `devtools' package. 

```{r}
install.packages("devtools")
library(devtools)
install_github("Muunraker/NIPALS-MCIA",force = TRUE)

```
The MATLAB version of NIPALS-MCIA is installed simply by downloading the files from Github 
and adding the `Functions' folder to the MATLAB path.

```{matlab}
addpath('<path to MATLAB_MCIA>\Functions')

```

## Using NIPALS-MCIA in R

The package currently includes only one vignette - a cut version of a three-omic NCI-60
cancer cell line dataset. Only one function call is needed to perform MCIA: 
```{r}
data(NCI60) # import data as "data_blocks"
mcia_results <- nipals_multiblock(data_blocks, preprocMethod='colprofile',num_PCs = 8, tol=1e-9)

```
Note: this dataset is reproduced from the omicade4 package (Meng et. al., 2014).

Meng, C., Kuster, B., Culhane, A. C., & Gholami, A. M. (2014). A multivariate approach to the integration of multi-omics datasets. BMC bioinformatics, 15(1), 1-13.

## Using NIPALS-MCIA in MATLAB

[Todo]
