# Function to mutate columns in row_data or col_data

Modification(mod) function to mutate data, row_data, col_data and
results. It is a wrapper function to dplyr::mutate() made for
metime_analyser objects.

## Usage

``` r
mod_mutate(object, which_data, type = "data", ...)
```

## Arguments

- object:

  An S4 object of class metime_analyser

- which_data:

  character or numeric input to define Dataset/result of interest. Has
  to be of length=1

- type:

  character input of length 1 to define the type of data to be
  manipulated. Accepted inputs are "row_data", "col_data", "data" and
  "results". However renamed results will be returned to the user as a
  list of results and will not return the full analyser object

- ...:

  arguments to pass directly into dplyr::mutate() function.

## Value

object with mutated data

## See also

[mod_filter](https://compneurobio.github.io/MeTime/reference/mod_filter.md),
[mod_rename](https://compneurobio.github.io/MeTime/reference/mod_rename.md)
