# Modification (mod) function to merge data and row_data (sample info) partially or completely.

a modification function that merges partial or complete data and
row_data of a dataset. If you want to use certain columns for the rest
of the analysis properly for annotating plots then use mod_mutate
function to rename the colnames and rownames (if applicable) to maintain
consistency in the merged dataset

## Usage

``` r
mod_merge_row_data_and_data(
  object,
  which_data,
  cols_list = list(data = NULL, row_data = NULL),
  append = FALSE,
  name = "merged_data"
)
```

## Arguments

- object:

  an S4 object of class metime_analyser.

- which_data:

  a character defining which data and row_data should be merged. Has to
  contain only one value.

- cols_list:

  A list of named character vectors to merge dataset. example:
  which_data = "dataset1"; cols_list = list(data=c(...),
  row_data=c(...)) or list(c(...), c(...)) If you want full data set
  cols_list=list(data=NULL, row_data=NULL) which is the default setting

- name:

  a character to define the name of a new dataset.

## Value

a new S4 object of class metime_analyser with the new merged dataset
appended to it

## See also

[mod_merge_results](https://compneurobio.github.io/MeTime/reference/mod_merge_results.md),
[mod_merge_data](https://compneurobio.github.io/MeTime/reference/mod_merge_data.md)
