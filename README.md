
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nipalsMCIA: Software to Compute Multi-Block Dimensionality Reduction

<!-- badges: start -->

[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/nipalsMCIA.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/nipalsMCIA)
[![R-CMD-check](https://github.com/Muunraker/nipalsMCIA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Muunraker/nipalsMCIA/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This package computes Multiple Co-Inertia Analysis (MCIA) on multi-block
data using the Nonlinear Iterative Partial Least Squares (NIPALS)
method.

Features include:

- Efficient computation of deflation and variance enabling embedding of
  high-volume (e.g. single-cell) datasets.  
- Functionality to perform out-of-sample embedding.
- Easy-to-adjust options for deflation and pre-processing
- Multiple visualization and analysis options for sample- and
  feature-level embedding results
- Streamlined and well-documented and supported code that is consistent
  with published theory to enable more efficient algorithm development
  and extension

**Citation**

For more information on the methodology used in nipalsMCIA and to cite,
please see

- Max Mattessich, Joaquin Reyna, Edel Aron, Ferhat Ay, Misha Kilmer, Steven H Kleinstein, Anna Konstorum, nipalsMCIA: flexible multi-block dimensionality reduction in R via nonlinear iterative partial least squares, Bioinformatics, Volume 41, Issue 1, January 2025, btaf015, https://doi.org/10.1093/bioinformatics/btaf015

## Installation

This package can be installed [via
Bioconductor](https://bioconductor.org/packages/release/bioc/html/nipalsMCIA.html):

``` r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("nipalsMCIA")
```

You can install the development version of nipalsMCIA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Muunraker/nipalsMCIA", ref = "devel",
                         force = TRUE, build_vignettes = TRUE)
```

## Basic Example

The package currently includes one test dataset: `data_blocks`. This is
a list of dataframes containing observations of variables from three
omics types (mRNA, proteins, and micro RNA) on 21 cancer cell lines from
the NCI60 cancer cell lines. The data file includes a `metadata` data
frame containing the cancer type associated with each cell line.

``` r
# load the package and set a seed for reproducibility
library(nipalsMCIA)
set.seed(42)
```

``` r
data(NCI60) # import data as "data_blocks" and metadata as "metadata_NCI60"

# examine the data and metadata
summary(data_blocks)
#>       Length Class      Mode
#> mrna  12895  data.frame list
#> miRNA   537  data.frame list
#> prot   7016  data.frame list
head(metadata_NCI60)
#>            cancerType
#> CNS.SF_268        CNS
#> CNS.SF_295        CNS
#> CNS.SF_539        CNS
#> CNS.SNB_19        CNS
#> CNS.SNB_75        CNS
#> CNS.U251          CNS
table(metadata_NCI60)
#> cancerType
#>      CNS Leukemia Melanoma 
#>        6        6        9
```

*Note: this dataset is reproduced from the [omicade4
package](https://www.bioconductor.org/packages/release/bioc/html/omicade4.html)
(Meng et. al., 2014). This package assumes all input datasets are in
sample by feature format.*

The main MCIA function can be called on `data_blocks` and optionally can
include `metadata_NCI60` for plot coloring by cancer type:

``` r
# to convert data_blocks into an MAE object we provide the simple_mae() function
data_blocks_mae <- simple_mae(data_blocks, row_format = "sample",
                              colData = metadata_NCI60)

mcia_results <- nipals_multiblock(data_blocks_mae = data_blocks_mae,
                                  col_preproc_method = "colprofile",
                                  num_PCs = 10, tol = 1e-12,
                                  color_col = "cancerType")
```

<img src="man/figures/README-call-mcia-1.png" width="100%" style="display: block; margin: auto;" />

Here `num_PCs` is the dimension of the low-dimensional embedding of the
data chosen by the user.
