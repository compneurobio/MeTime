# Meta comparison for regression outputs

Compare regression outputs within or across results.

## Usage

``` r
meta_regression(
  object,
  method = c("sign", "cor", "het"),
  result_index = NULL,
  name = "meta_regression_1"
)
```

## Arguments

- object:

  a S4 object of class metime_analyser or a list of two metime_analyser
  objects

- method:

  a character vector of methods 'sign', 'cor', 'het'

- result_index:

  character/numeric input for results. If NULL, all matching results are
  used.

- name:

  a character input to set the name of the results

## Value

An S4 object of class meta_results with the compared results and meta
results
