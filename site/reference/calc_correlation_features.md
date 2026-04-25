# Function to calculate correlation

calculate pairwise correlations This function creates a dataframe for
plotting from a dataset.

## Usage

``` r
calc_correlation_features(
  object,
  which_data,
  method,
  name = "calc_correlation_pairwise_1",
  stratifications
)
```

## Arguments

- object:

  S4 Object of class metime_analyser

- which_data:

  specify datasets to calculate on. One or more possible

- method:

  default setting: method="pearson", Alternative "spearman" also
  possible

- name:

  name of the results should be of length=1

- stratifications:

  List to stratify data into a subset. Usage list(name=value)

## Value

data.frame with pairwise results

## Examples

``` r
# Example to calculate correlations
dist <- calc_correlation(object=metime_analyser_object, which_data="name of the dataset", 
          method="pearson")
#> Error in calc_correlation(object = metime_analyser_object, which_data = "name of the dataset",     method = "pearson"): could not find function "calc_correlation"
```
