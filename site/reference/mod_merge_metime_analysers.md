# Function to merge one or more metime_analyser objects

function to merge multiple metime_analyser objects. Will not be
displayed in add_function_info()

## Usage

``` r
mod_merge_metime_analysers(..., annotations_index)
```

## Arguments

- ...:

  metime analyser objects that are to be merged. The objects can be
  named by parsing it as an argument If not named then they will be
  named by the function

- annotations_index:

  new list with annotations_index. If set to NULL the first object's
  annotations are taken with modification.

## Value

A merged metime_analyser object

## See also

[add_function_info](https://compneurobio.github.io/MeTime/reference/add_function_info.md)
