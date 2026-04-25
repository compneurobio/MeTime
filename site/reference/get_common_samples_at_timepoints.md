# Get common samples at multiple timepoints chosen

Get a vector of common samples that can be used for stratification.

## Usage

``` r
get_common_samples_at_timepoints(object, which_data, timepoints)
```

## Arguments

- which_data:

  a character to define which dataset is to be used.

- timepoints:

  a character vector to define time points

- objecta:

  S4 object of the class "metime_analyzer".

## Value

A character vector of common subjects at common time points.
