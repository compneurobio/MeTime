# Summarize datasets in a metime_analyser object

Returns dataset size and missingness summaries.

## Usage

``` r
summarize_dataset(object, which_data = NULL)
```

## Arguments

- object:

  a S4 object of class "metime_analyser".

- which_data:

  character vector of dataset names to summarize. If NULL, summarizes
  all datasets.

## Value

A data.frame with summary statistics for each dataset.
