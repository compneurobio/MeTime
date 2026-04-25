# Function to calculate dissimilarity using distance measures

calculate pairwise distances between samples This function creates a
dataframe for plotting from a dataset.

## Usage

``` r
calc_distance_samples(
  object,
  which_data,
  method,
  name = "calc_distance_pairwise_1",
  stratifications
)
```

## Arguments

- object:

  S4 Object of class metime_analyser

- which_data:

  specify datasets to calculate on. One or more possible

- method:

  default setting: method="euclidean", Alternative "maximum","minimum",
  "manhattan","canberra","minkowski" are also possible

- name:

  name of the results should be of length=1

- stratifications:

  List to stratify data into a subset. Usage list(name=value)

## Value

data.frame with pairwise results

## Examples

``` r
# Example to calculate pairwise distances
dist <- calc_pairwise_distance(object=metime_analyser_object, which_data="name of the dataset", 
          method="euclidean")
#> Error in calc_pairwise_distance(object = metime_analyser_object, which_data = "name of the dataset",     method = "euclidean"): could not find function "calc_pairwise_distance"
```
