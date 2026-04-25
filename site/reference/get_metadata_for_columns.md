# Get metadata for columns(in most cases for metabolites)

function to generate a metadata list for building the MeTime plotter
object

## Usage

``` r
get_metadata_for_columns(object, which_data, columns = NULL)
```

## Arguments

- object:

  a S4 object of the class "metime_analyzer".

- which_data:

  a character to define which dataset is to be used.

- columns:

  A named list of a vector of named characters containing the columns of
  interest. Length of the list should be same as length of which_data.
  List can be named or unnamed but the character vectors should be
  named. Moreover, they should have same names. Default is set to NULL
  which results in empty metadata dataframe. Ex:
  list(nmr_data=c(id="id", sub_pathway="Group"), lipid_data=c(id="id",
  sub_pathway="sub_pathway"))

## Value

a dataframe with metadata information.
