# Function to calculate fold change

Function to calculate relative log transformed concentration of
metabolites with respect to baseline concentration (foldchange) of
multiple datasets at once

## Usage

``` r
mod_trans_foldchange(object, which_data, log2 = T, col = "all")
```

## Arguments

- object:

  an S4 object of class metime_analyser.

- which_data:

  a character vector defining the name of the dataset/s to be used.

- log2:

  logical input. TRUE denotes that the data is already log transformed.

- col:

  character vector to select the metabolites. Set to "all" for all
  metabolites.

## Value

a S4 object of the class metime_analyzer with datasets updated with
foldchange values.
