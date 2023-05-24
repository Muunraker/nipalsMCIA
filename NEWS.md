# 0.99.3 (2023-05-26):

## Major changes

* Changed eigenvalue calculation in `nipals_iter()` to variance of global scores instead of directly computing the singular value.

## Minor improvements and bug fixes

* Changed vignette styling from `rmdformats` to `BiocStyle`
* Removed empty helper.R
* Added significantly more unit testing
* Fixed bug in `projection_plots()` when metadata was provided but no `color_col` was selected.
* Added checks for consistency in sample names across data blocks and metadata. 

# 0.99.2 (2023-03-25):

## Major changes

* Shrank the vignettes sizes (especially Vignette 2)

## Minor improvements and bug fixes

* Restructured Vignette 2 to be more streamlined and have more explanations
* Add an additional single cell data file to the repository using `piggyback`

# 0.99.1 (2023-02-26):

## Major changes

* Included support for `MultiAssayExperiment` in `nipals_multiblock`
* Improved access to existing data objects

## Minor improvements and bug fixes

* Made `get_colors()` more flexible for different color palette options

# 0.99.0 (2022-10-21):

* Added single cell data to the repository using `piggyback`
