# Convert SummarizedExperiment objects back to metime_analyser

This helper reverses the behaviour of `get_convert_to_se` by creating a
new `metime_analyser` from one or more `SummarizedExperiment` objects.
It reconstructs the original data matrix orientation (samples as rows
and features as columns) and retrieves row- and column-level metadata
stored in the `colData` and `rowData` of the SummarizedExperiment,
respectively. Datasets that lack the mandatory `id`, `subject` and
`time` columns in their row-level metadata will trigger an error, as
these fields are required for a valid `metime_analyser`.

## Usage

``` r
get_convert_from_se(se, name = NULL, annotations_index = list())

# S4 method for class 'list'
get_convert_from_se(se, name = NULL, annotations_index = list())

# S4 method for class 'MultiAssayExperiment'
get_convert_from_se(se, name = NULL, annotations_index = list())
```

## Arguments

- se:

  A single `SummarizedExperiment` or a list of such objects. When a list
  is provided, each entry will be added as a separate dataset to the
  returned `metime_analyser` object.

- name:

  Optional character vector of dataset names. If `se` is a list and
  `name` is `NULL`, the names of the list elements will be used. When
  `se` is a single object and `name` is `NULL`, the default dataset name
  "set_1" will be applied.

- annotations_index:

  A named list specifying phenotype and medication datasets. This
  mirrors the argument in `get_make_analyser_object`. When converting
  multiple datasets, the same annotations list will be applied to all.

## Value

An object of class `metime_analyser` containing one or more
reconstructed datasets.

## Functions

- `get_convert_from_se(list)`: Convert a list of SummarizedExperiment
  objects. When a single dataset has been split by timepoint (i.e. the
  list contains multiple SummarizedExperiments whose names share a
  common prefix before the final underscore), this method automatically
  combines them back into a single dataset by row-binding their
  sample-level data and taking the feature-level metadata from the first
  element.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Assume 'se_obj' is a SummarizedExperiment produced by get_convert_to_se()
  analyser <- get_convert_from_se(se_obj, name = "my_dataset")
} # }
```
