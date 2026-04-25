# Filter samples by missingness

Removes samples with missingness above a threshold.

## Usage

``` r
mod_filter_samples_by_missingness(object, threshold = 0.2, which_data = NULL)
```

## Arguments

- object:

  a S4 object of class "metime_analyser".

- threshold:

  numeric between 0 and 1 defining maximum allowed missingness.

- which_data:

  character vector of dataset names to filter. If NULL, filters all
  datasets.

## Value

A modified metime_analyser object.
