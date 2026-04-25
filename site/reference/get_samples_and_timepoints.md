# Summarize number of time points and the total number of samples available at that point

Obtain the number of unique samples available at each time point.

## Usage

``` r
get_samples_and_timepoints(object, which_data)
```

## Arguments

- which_data:

  a character to define which dataset is to be used.

- objecta:

  a S4 object of the class "metime_analyzer".

## Value

A dataframe with time points and number of samples at each time point.
