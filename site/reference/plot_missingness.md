# Plot missingness heatmap for a dataset

Creates a heatmap of missing values for a dataset.

## Usage

``` r
plot_missingness(object, which_data, max_features = 100)
```

## Arguments

- object:

  a S4 object of class "metime_analyser".

- which_data:

  a single dataset name.

- max_features:

  maximum number of features to include in the plot.

## Value

A ggplot object.
