# 1.3.2 (2024-12-10)

## Major changes

## Minor improvements and bugfixes
* Added `harmonize = FALSE` option to `nipals_multiblock()` to allow disabling (slow) sample harmonization for large datasets.
* Fixed a bug in `nipals_multiblock()` to allow the `numPCs =1` argument (i.e. ability to compute a single global score).

# 1.2.1 (2024-08-31)

## Major changes

## Minor improvements and bugfixes
* Updated 'Single-Cell-Analysis' vignette to remove potential issue with BiocFileCache having two copies of data. 
* Updated author list and citation. 
* Updated references in readme. 

# 1.2.0 (2024-04-26)

## Major changes

## Minor improvements and bugfixes
* Clarified the dataset used in the single-cell vignette and reformatted the workflow diagram
* Removed the final log transformation in the single-cell vignette when saving for MCIA
* Minor plotting improvements (removed redundancies in legend names, consistency of axis text, etc.)

---

# 1.1.0 (2024-03-21)

## Major changes
* Added `BiocFileCache` to replace previous method of downloading large SC datasets. 

## Minor improvements and bugfixes
* Updated `nipals_iter` to only use `var(gs)` to compute the significance of global scores.
* Added to documentation for `predict_gs()` warning against use with CPCA deflation. 

---

# 0.99.7 (2023-10-07):

## Minor improvements and bug fixes

* Updated README to reflect MAE changes.
* Updated the citation. 
* Fixed bug in documentation for `nipals_multiblock`.
* Stopped downloading one of the files for the single cell analysis vignette.

---

# 0.99.6 (2023-09-09):

## Major changes

* Switched primary input to `nipals_multiblock()` to a `MultiAssayExperiment` object.
* `nipals_multiblock()` now outputs an object of the `NipalsResult` class.
* Converted all downstream analysis functions to work with the `NipalsResult` class.

## Minor improvements and bug fixes

* Added `simple_mae()` function to convert a list of dataframes to a `MultiAssayExperiment` object.
* Fixed missing `\value` fields in man page for `NipalsResult` class. 

---

# 0.99.5 (2023-06-22):

## Major changes

* Bumping version number to trigger re-build following bug in BiocCheck.

---

# 0.99.4 (2023-06-01):

## Major changes

* Fixed data corruption issues from v0.99.3.

---

# 0.99.3 (2023-06-01):

## Major changes

* Changed the eigenvalue calculation in `nipals_iter()` to compute the variance of the global score at each deflation step. Prior versions used an SVD method to compute the singular values of the deflated data matrix directly.  
* Fixed bug in `projection_plot()` where there was a mismatch between color labels and plotting order.
* Added parameter to `nipals_multiblock` that specifies whether the samples are in the rows or the columns.

## Minor improvements and bug fixes

* Changed vignette styling from `rmdformats` to `BiocStyle` and added installation sections to all of the vignettes.
* Removed empty helper.R file.
* Added significantly more unit testing.
* Fixed bug in `projection_plot()` when metadata was provided but no `color_col` was selected.
* Renamed the associated output of `col_preproc_method` in `nipals_multiblock`. The metadata field is also now available in the output independent of whether metadata is provided in the input.
* Added checks for consistency in sample names across data blocks and metadata. 

---

# 0.99.2 (2023-03-25):

## Major changes

* Shrank the vignettes sizes (especially the single cell analysis vignette).

## Minor improvements and bug fixes

* Restructured the single cell analysis vignette to be more streamlined and have more explanations.
* Add an additional single cell data file to the repository using `piggyback`.

---

# 0.99.1 (2023-02-26):

## Major changes

* Included support for `MultiAssayExperiment` in `nipals_multiblock`.
* Improved access to existing data objects.

## Minor improvements and bug fixes

* Made `get_colors()` more flexible for different color palette options.

---

# 0.99.0 (2022-10-21):

* Added single cell data to the repository using `piggyback`.
