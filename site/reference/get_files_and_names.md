# Pack a dataset (data, row_data, col_data) into a single object of class "metime_analyser".

This function loads all the files from the parent directory. It assumes
a certain naming pattern as follows: "datatype_None\|col\|row_data.rds"
Any other naming pattern is not allowed. The function first writes all
files into a list and each type of data is packed into its respective
class i.e. col_data, row_data or data.

## Usage

``` r
get_files_and_names(path, annotations_index)
```

## Arguments

- path:

  a character defining the path to the parent directory

- annotations_index:

  a named list to be filled as list(phenotype="Name or index of the
  files", medication="Name or index of the files")

## Value

An object of class "metime_analyser".

## Examples

``` r
get_files_and_names(path="/path/to/parent/directory", 
  annotations_index=list(phenotype="Name of phenotype file", 
medication="name of phenotype file"))
#> Datasets: 
#> Results: 
```
