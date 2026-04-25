# Modification (mod) function to merge two sets of data or partial data merging for any analysis

a modification function that merges two sets of data or partial data
including data, col_data and row_data. If you want to use certain
columns for the rest of the analysis properly for annotating plots then
use mod_mutate function to rename the colnames and rownames (if
applicable) to maintain consistency in the merged dataset

## Usage

``` r
mod_merge_data(
  object,
  which_data,
  filter_samples = NULL,
  cols_list = NULL,
  append = FALSE,
  name = "merged_data"
)
```

## Arguments

- object:

  an S4 object of class metime_analyser.

- which_data:

  a vector of character defining which data should be merged. Has to
  contain two or more values. If you want to merge partial data from
  other datasets then the first dataset will be consided to be complete
  unless and until you specify the first dataset as name in the
  cols_for_meta list. See below for example

- filter_samples:

  a character specifying if samples should be filtered. Default set to
  no filtering. Other options include common samples (samples used if
  found in all datasets), or set to name of dataset (samples filtered by
  samples of one dataset)

- cols_list:

  A list of named character vectors to merge dataset. If the first
  dataset's name is used in the list then those columns will be
  subsetted else the order of which_data will be followed irrespective
  of the names example: which_data = c("dataset1", "dataset2");
  cols_list = list(dataset1=c(...), dataset2=c(...)) or list(c(...)) If
  you want full data set cols_list=NULL

- append:

  logical. TRUE adds this new dataset to the object. FALSE creates a new
  object with just the data and the latest results section. Default is
  set to FALSE

- name:

  a character to define the name of a new dataset.

## Value

a new S4 object of class metime_analyser with the new merged dataset
appended to it

## See also

[mod_merge_results](https://compneurobio.github.io/MeTime/reference/mod_merge_results.md),
[mod_merge_row_data_and_data](https://compneurobio.github.io/MeTime/reference/mod_merge_row_data_and_data.md)
