# Imputation of missing values

This function imputes missing values using one of the implemented
imputation methods.

## Usage

``` r
mod_impute(
  object,
  which_data,
  method = "rf",
  thresh = 0.3,
  path_to_diagnostics = NULL
)
```

## Arguments

- object:

  a S4 object of class "metime_analyser".

- which_data:

  a character defining which dataset should be imputed.

- method:

  a character defining the method for imputation. See details for more
  information. Default set to 'rf'

- thresh:

  a numeric defining at which missingness level metabolites will be
  excluded. Default set to 0.3 (30% missingness)

- path_to_diagnostics:

  a character defining the path to where a diagnostic file is being
  written.

## Value

S4 object of the class "metime_analyser" with mutated col_data and
row_data.

## Details

additional methods are currently not implemented such as knn, pmm, norm
etc.
