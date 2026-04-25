# This function appends an object of class metime_analyser with a new dataset.

function to apply on metime_analyser object to append a new dataset into
the existing object

## Usage

``` r
add_dataset(object, data, col_data, row_data, name)
```

## Arguments

- object:

  S4 object of class metime_analyser

- data:

  data.frame containing data

- col_data:

  data.frame containing col_data: id column of col data has to match
  colnames of data

- row_data:

  data.frame containing row_data: id column of row data has to match
  rownames of data

- name:

  Name of the new dataset

## Value

An object of class metime_analyser

## Examples

``` r
# append data frames into the metime_analyser object
appended_object <- add_dataset(object=metime_analyser_object, 
                   data=data, row_data=row_data, col_data=col_data, name="name of the new dataset")
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'add_dataset': object 'metime_analyser_object' not found
```
