# Function to update plots post calculations

Modification(mod) Function to generate all possible/available plots of a
calculation. This function is a wrapper function for plot()

## Usage

``` r
mod_generate_plots(object, .interactive = FALSE, type, results_index = NULL)
```

## Arguments

- object:

  An S4 object of class metime_analyser

- .interactive:

  logical to make the plot interactive or not

- type:

  character to define the type of calculation used for updating the
  plot. Allowed arguments are "network", "PCA", "UMAP", "tSNE",
  "CI_metabotype", "CI_metabolite", "regression", "pairwise_distance",
  "pairwise_correlation", "colinearity", "regression",
  "distribution_samples, "distribution_metabs", "feature_selection",
  "clusters" cols_for_samples, cols_for_metabs, cols_for_meta etc will
  be used. So make sure you set those correctly for better results.
  Allowed inputs are: c("ggm\|network", "PCA", "UMAP", "tSNE",
  "CI_metabotype", "CI_metabolite", "pairwise_distance",
  "pairwise_correlation", "colinearity", "distribution_metabs",
  "distribution_samples", "regression", "feature_selection")

- results_index:

  character/numeric input to define the results that you want to plot or
  replot with our automation. Default will be set to NULL. Length of
  results_index should be equal to 1.

## Value

object with plots of the newest calculation

## See also

[plot](https://rdrr.io/r/graphics/plot.default.html)
