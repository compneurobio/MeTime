# Function to calculate colinearity of features

Function to calculate colinearity in a dataset

## Usage

``` r
calc_colinearity_features(
  object,
  which_data,
  name = "calc_colinearity_1",
  stratifications = list()
)
```

## Arguments

- object:

  An S4 object of class metime_analyser

- which_data:

  Dataset to check for colinearity

- name:

  character to define the name of the result

- stratifications:

  List to stratify data into a subset. Usage list(name=values)

- show_all:

  logical. True will only filter out colinear data

## Value

plotter object with data for heatmap information
