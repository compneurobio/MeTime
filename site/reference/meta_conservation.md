# Meta comparison for conservation index results

Compare conservation index results within a result or across results.

## Usage

``` r
meta_conservation(
  object,
  result_index = NULL,
  top_k = 50,
  name = "meta_conservation_1"
)
```

## Arguments

- object:

  a S4 object of class metime_analyser or a list of two metime_analyser
  objects

- result_index:

  character/numeric input for results. If NULL, all matching results are
  used.

- top_k:

  numeric indicating the top-K features used for overlap calculations

- name:

  a character input to set the name of the results

## Value

An S4 object of class meta_results with the compared results and meta
results
