# Calculate significance threshold after correcting for independent tests.

Eigenvalue based calculation of independent test (removing colinearity)
described by Li et al. (See details for more)

## Usage

``` r
get_li_thresh(object, which_data, verbose = F)
```

## Arguments

- object:

  a S4 object of the class "metime_analyzer".

- which_data:

  a character to define which dataset is to be used.

- verbose:

  a logical to print out the number of independent tests. Default is set
  to FALSE.

## Value

a numeric value describing the Li threshold.

## Details

A detailed description of the method can be found in [Li et al. 2012,
Evaluating the effective numbers of independent tests and significant
p-value thresholds in commercial genotyping arrays and public imputation
reference datasets](https://doi.org/10.1007%2Fs00439-011-1118-2)
