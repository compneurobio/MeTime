# Get results list for S4 object of class "metime_analyser" object

Compile a result element for metime_analyser object.

## Usage

``` r
get_make_results(object, data, metadata, calc_type, calc_info, name)
```

## Arguments

- object:

  a S4 object of the class "metime_analyzer".

- data:

  a list of dataframes of plotable data obtained from any calc function.

- metadata:

  a dataframe or a list of dataframes with the metadata for the plot
  table mentioned above. To obtain these see get_metadata_for_rows() and
  get_metadata_for_columns().

- calc_type:

  a character vector to specify type of calculation. Should be the same
  length as the list data provided.

- calc_info:

  a character to define the information about calculation, should be the
  same length as the list data provided.

- name:

  a character to be used as the index of the result.

## Value

A S4 object of class "metime_analyser" with the result appended to the
list of results.
