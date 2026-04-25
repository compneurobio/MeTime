# Add a data frame to analyzer object.

Add a data.frame x to existing data within an analyzer object. Data can
be added to data, col_data or row_data. A full join of the data is used.

## Usage

``` r
add_data(object, which_data, type, x, id)
```

## Arguments

- object:

  a S4 object of the class "metime_analyzer".

- which_data:

  a character to define which dataset is to be used.

- type:

  a character defining the data which should be appended. Can only be
  data, col_data, row_data. Default set to data.

- x:

  a data.frame that is merged to the data.

- id:

  a character vector defining the column to be used for matching.
  Default set to NULL (first column of x will be used as id column)

## Value

a S4 object of class "metime_analyser" with appended data.frame x
