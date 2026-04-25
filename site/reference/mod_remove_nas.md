# Function to remove NA's from data matrices

Modification(mod) function applied on S4 object to remove NA's and
change col_data and row_data accordingly

## Usage

``` r
mod_remove_nas(object, which_data)
```

## Arguments

- object:

  S4 object of class metime_analyser

- which_data:

  character vector to define the datasets to be used. Can also use a
  numeric vector with indices.

## Value

S4 object with NA's removed and col_data, row_data manipulated
accordingly
