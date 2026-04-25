# Modification function to select columns in data, row_data or col_data

Modification (mod) function to select columns in data, row_data or
col_data. It is a wrapper function to dplyr::select() made for
metime_analyser objects.

## Usage

``` r
mod_select(object, which_data, type = "data", ...)
```

## Arguments

- object:

  A S4 object of class metime_analyser.

- which_data:

  a vector of character defining which data/result is to be filtered.

- type:

  character input of length 1 to define the type of data to be
  manipulated. Accepted inputs are "row_data", "col_data", "data" and
  "results". However renamed results will be returned to the user as a
  list of results and will not return the full analyser object

- ...:

  arguments to pass directly into dplyr::filter() function.

## Value

object with mutated data, col_data or row_data. However, if type is
"results" then plot_data with modifications expected will be returned

## See also

[mod_mutate](https://compneurobio.github.io/MeTime/reference/mod_mutate.md),
[mod_rename](https://compneurobio.github.io/MeTime/reference/mod_rename.md),
[mod_filter](https://compneurobio.github.io/MeTime/reference/mod_filter.md)
