# Calculate baseline change per subject

Computes baseline-adjusted values for each sample.

## Usage

``` r
calc_baseline_change(
  object,
  which_data,
  baseline_timepoint,
  method = "diff",
  name = NULL
)
```

## Arguments

- object:

  a S4 object of class "metime_analyser".

- which_data:

  dataset name to use.

- baseline_timepoint:

  timepoint used as baseline.

- method:

  "diff" for subtraction or "ratio" for division.

- name:

  name for the results entry.

## Value

A metime_analyser object with results appended.
