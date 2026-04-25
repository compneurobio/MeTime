# Function for Plotting distributions of phenotypic variables

A method to be applied onto s4 object so as to obtain distributions of
various phenotypic variables

## Usage

``` r
calc_distribution_metabs(
  object,
  which_data,
  cols,
  name = "calc_distribution_metabs_1"
)
```

## Arguments

- object:

  An object of class metime_analyser

- which_data:

  Name of the dataset from which the samples will be extracted

- cols:

  character vector to define the columns whose distributions are wanted
  from col_data

- name:

  Character to name the results

## Value

S4 object with updated plot_data and plots with plots being either 1)
density plot 2) bar plot and a line plot
