# Add values measured at a screening time for samples to be added to all time points

Add all data points that were measured only during a defined screening
time point to all other measurements of the same subject. Note that this
function is only valid in the case of phenotype data and not in the case
of row_data

## Usage

``` r
add_screening_vars(object, vars)
```

## Arguments

- object:

  a S4 object of class "metime_analyser".

- vars:

  a character vector naming the columns of interest

## Value

S4 object of class "metime_analyser" with screening measurements added
to all other time points of the same subject.
