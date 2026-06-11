# Convert analyser datasets to SummarizedExperiment objects with optional time splitting

Converts each dataset in a `metime_analyser` into one or more
`SummarizedExperiment` objects. When `split_by_time` is `TRUE`,
longitudinal datasets (those whose row-level metadata include a `time`
column with more than one unique value) are split into separate
`SummarizedExperiment`s for each timepoint. Otherwise, the entire
dataset is stored in a single `SummarizedExperiment` with the `time`
variable preserved in the `colData` slot. Feature metadata (`rowData`)
is not duplicated across timepoints; the same `DataFrame` is reused for
each split.

Phenotype and medication datasets are skipped by default.

## Usage

``` r
get_convert_to_se(
  object,
  which_data = NULL,
  exclude = c("phenotype_data", "medication_data"),
  split_by_time = TRUE
)
```

## Arguments

- object:

  A `metime_analyser` object.

- which_data:

  Character vector of dataset names to convert. If `NULL`, all datasets
  are converted.

- exclude:

  Character vector of dataset names to skip. Defaults to
  `c("phenotype_data", "medication_data")`.

- split_by_time:

  Logical indicating whether longitudinal datasets should be split by
  timepoint. If `FALSE`, a single `SummarizedExperiment` is returned per
  dataset with a `time` column in the sample-level metadata.

## Value

A named list of `SummarizedExperiment` objects.
