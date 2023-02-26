
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nipalsMCIA: Software to Compute Multi-Block Dimensionality Reduction

<!-- badges: start -->

<!-- badges: end -->

This package computes Multiple Co-Inertia Analysis (MCIA) on multi-block
data using the Nonlinear Iterative Partial Least Squares (NIPALS)
method.

Features include:

  - Efficient computation of deflation and variance enabling embedding
    of high-volume (e.g. single-cell) datasets.  
  - Functionality to perform out-of-sample embedding.
  - Easy-to-adjust options for deflation and pre-processing
  - Multiple visualization and analysis options for sample- and
    feature-level embedding results
  - Streamlined and well-documented and supported code that is
    consistent with published theory to enable more efficient algorithm
    development and extension

**References**

Mattesich (2022) A Review of Multi-Block Dimensionality Reduction via
Multiple Co-Inertia Analysis, M.S. Thesis, Dept. of Mathematics, Tufts
University (<https://dl.tufts.edu/concern/pdfs/cz30q6773>)

Hanafi et al. (2011) Connections between multiple co-inertia analysis
and consensus principal component analysis, Chemometrics and Intelligent
Laboratory Systems 106 (1)
(<https://doi.org/10.1016/j.chemolab.2010.05.010>.)

Meng et al. (2014) A multivariate approach to the integration of
multi-omics datasets, BMC Bioinformatics 2014(15)
(<https://doi.org/10.1186/1471-2105-15-162>)

## Installation

This package currently can only be installed using
`devtools::install_github()`. A CRAN/Bioconductor version is in
progress.

You can install the development version of nipalsMCIA from
[GitHub](https://github.com/) with:

``` r
library(devtools)
devtools::install_github("Muunraker/nipalsMCIA", ref = "code-development",
                         force = TRUE, build_vignettes = TRUE)
```

``` r
library(nipalsMCIA)
```

## Basic Example

The package currently includes one test dataset: `data_blocks`. This is
a list of dataframes containing observations of variables from three
omics types (mRNA, proteins, and micro RNA) on 21 cancer cell lines from
the NCI60 cancer cell lines. The data file includes a `metadata` data
frame containing the cancer type associated with each cell line.

``` r
data(NCI60) # import data as "data_blocks" and metadata as "metadata_NCI60"

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
#> metadata_NCI60
#>      CNS Leukemia Melanoma 
#>        6        6        9
```

Note: this dataset is reproduced from the [omicade4
package](https://www.bioconductor.org/packages/release/bioc/html/omicade4.html)
(Meng et. al., 2014). This package assumes all input datasets are in
sample by feature format.

The main MCIA function can be called on `data_blocks` and optionally can
include `metadata_NCI60` for plot coloring by cancer type:

``` r
set.seed(42)
mcia_results <- nipals_multiblock(data_blocks, preproc_method = 'colprofile',
                                  metadata = metadata_NCI60, color_col = "cancerType", 
                                  num_PCs = 10, tol = 1e-12)
```

<img src="man/figures/README-call-mcia-1.png" width="100%" style="display: block; margin: auto;" />

Here `numPCs` is the dimension of the low-dimensional embedding of the
data chosen by the user.
