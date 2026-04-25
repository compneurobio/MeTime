# Filter timepoints from datasets

Removes samples that are not in the specified timepoints.

## Usage

``` r
mod_filter_timepoints(object, timepoints, which_data = NULL)
```

## Arguments

- object:

  a S4 object of class "metime_analyser".

- timepoints:

  character or numeric vector of timepoints to keep.

- which_data:

  character vector of dataset names to filter. If NULL, filters all
  datasets.

## Value

A modified metime_analyser object.
