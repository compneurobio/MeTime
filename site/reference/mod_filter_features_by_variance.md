# Filter features by variance

Removes features with variance below a threshold.

## Usage

``` r
mod_filter_features_by_variance(object, min_variance = 0, which_data = NULL)
```

## Arguments

- object:

  a S4 object of class "metime_analyser".

- min_variance:

  numeric defining minimum variance required to keep a feature.

- which_data:

  character vector of dataset names to filter. If NULL, filters all
  datasets.

## Value

A modified metime_analyser object.
