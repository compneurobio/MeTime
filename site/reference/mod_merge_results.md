# Function to merge two sets of results

Modification(mod) function to merge two sets of results based on
calc_type

## Usage

``` r
mod_merge_results(object, results_index, sub_results = 1, groups, name)
```

## Arguments

- object:

  An S4 object of class metime_analyser

- results_index:

  character of numeric vector to define the index of results to merge.
  if length of results_index is 1, then the plot_data list is merged
  into a single data.frame

- sub_results:

  List of character/numeric vectors to define indices of plot_data to be
  merged. This is not needed when trying to merge plot_data of same
  results but need it for merging different results. The length of
  sub_results should be equal to length of results_index and also should
  be in the same order.

- groups:

  character vector to name the results you want to merge and should be
  of length plot_data in the case of merging results of same
  calculation. Else groups should be of length equal to
  unlist(sub_results) and the groups should follow the order with
  precedence to results_index. Eg: results_index=1:2,
  sub_results=list(c(3,4), c(5,6)), then groups \<- c("name_1_3",
  "name_1_4", "name_2_5", "name_2_6")

- name:

  character to name the new merged results

## Value

metime_analyser object with merged results

## See also

[mod_merge_data](https://compneurobio.github.io/MeTime/reference/mod_merge_data.md),
[mod_merge_row_data_and_data](https://compneurobio.github.io/MeTime/reference/mod_merge_row_data_and_data.md)
