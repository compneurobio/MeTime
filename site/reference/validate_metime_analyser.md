# Validate a metime_analyser object

Check consistency of list_of_data, row_data, and col_data.

## Usage

``` r
validate_metime_analyser(object, which_data = NULL)
```

## Arguments

- object:

  a S4 object of class "metime_analyser".

- which_data:

  character vector of dataset names to validate. If NULL, validates all
  datasets.

## Value

A list with dataset-level issues and a summary table.
