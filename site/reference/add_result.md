# Add a result element to a metime_analyser

Add a result from a calc\_ function to the metime_analyser.

## Usage

``` r
add_result(
  object,
  x,
  name = NULL,
  functions_applied = NULL,
  calc_info = NULL,
  calc_type = NULL,
  results_index = NULL
)
```

## Arguments

- object:

  a S4 object of class "metime_analyser"

- x:

  a data.frame or plot that should be added to the list of results

- functions_applied:

  a character providing information on the calculation

- calc_info:

  a character that provides information on the calculation

- calc_type:

  a character that defines the type of calculation

- results_index:

  a haracter or numeric of length=1 to define the results to which you
  want to add x. Default is set to NULL and that implies that new
  results are generated

## Value

S4 object of class "metime_analyser" with appended results
