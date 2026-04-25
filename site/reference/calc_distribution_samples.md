# Function for Plotting distributions of phenotypic variables

A method to be applied onto s4 object so as to obtain distributions of
various phenotypic variables

## Usage

``` r
calc_distribution_samples(
  object,
  which_data,
  cols,
  stratifications,
  name = "calc_distribution_samples_1"
)
```

## Arguments

- object:

  An object of class metime_analyser

- which_data:

  Name of the dataset from which the samples will be extracted

- cols:

  character vector to define the columns whose distributions are wanted
  from row_data

- stratifications:

  List to define the stratification of interest

- name:

  Character to name the results

## Value

S4 object with updated plot_data and plots with plots being either 1)
density plot 2) bar plot
