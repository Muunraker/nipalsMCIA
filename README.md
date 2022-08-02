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

[Todo] 


## Using NIPALS-MCIA in MATLAB

[Todo]
