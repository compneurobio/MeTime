# Get metadata for rows(in most cases for samples)

Generate a list of metadata for the S4 object of the class
"metime_analyzer".

## Usage

``` r
get_metadata_for_rows(object, which_data, columns)
```

## Arguments

- object:

  a S4 object of the class "metime_analyzer".

- which_data:

  a character to define which dataset is to be used.

- columns:

  A list of character vectors for the columns of interest. Length of the
  list should be same as length of which_data.

## Value

A data frame with metadata information for rows.
