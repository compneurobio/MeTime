# Function to rename columns in the datasets and results of analyser object

Function to rename columns of data, row_data, col_data and plot_data of
the analyser object. Arguments will be passed directly into
dplyr::rename() function

## Usage

``` r
mod_rename(object, which_data, type = "data", ...)
```

## Arguments

- object:

  An S4 object of class metime_analyser

- which_data:

  character or numeric input of length 1 to define the index of dataset
  or results

- type:

  character input of length 1 to define the type of data to be
  manipulated. Accepted inputs are "row_data", "col_data", "data" and
  "results". However renamed results will be returned to the user as a
  list of results and will not return the full analyser object

- ...:

  arguments to pass to dplyr::rename to change the name of the
  functions. Example new_colname = "old_colname"

## Value

metime_analyser object with mutated names of dataframes that are parsed

## See also

[mod_filter](https://compneurobio.github.io/MeTime/reference/mod_filter.md),
[mod_mutate](https://compneurobio.github.io/MeTime/reference/mod_mutate.md),
[mod_select](https://compneurobio.github.io/MeTime/reference/mod_select.md)
