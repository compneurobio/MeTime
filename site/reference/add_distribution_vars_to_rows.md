# Add phenotypic measurements that are not added to the row_data of the dataset

A method applied on the s4 object of class "metime_analyser" to add all
those datapoints that are present in phenotype dataframe but not in
row_data comes with the feature of updating those data points measured
only at screening to all datapoints and then adding it to row_data

## Usage

``` r
add_distribution_vars_to_rows(
  object,
  screening_vars,
  distribution_vars,
  which_data
)
```

## Arguments

- object:

  An object of class metime_analyser

- screening_vars:

  A character vector to define the vars that are to be updated as per
  add_screening_vars() else set it to NULL(Default).

- distribution_vars:

  A character naming the vars of interest

- which_data:

  dataset to which the information is to be added(only 1 can be used at
  a time)

## Value

object of class metime_analyser with phenotype data added to row data

## See also

[add_screening_vars](https://compneurobio.github.io/MeTime/reference/add_screening_vars.md)

## Examples

``` r
object <- add_distribution_vars_to_rows(object=data, screening_vars=c("var1", "var2"), 
    distribution_vars=c("var1", "var2", "var3"), which_data="dataset1")
#> Error: unable to find an inherited method for function ‘add_distribution_vars_to_rows’ for signature ‘object = "function"’
```
