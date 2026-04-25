# Convert analyser datasets to SummarizedExperiment objects

Converts each dataset in a metime_analyser into a SummarizedExperiment
object. Phenotype and medication datasets are skipped. Longitudinal
datasets are split into separate SummarizedExperiment objects by
timepoint with a timepoint suffix.

## Usage

``` r
get_convert_to_se(
  object,
  which_data = NULL,
  exclude = c("phenotype_data", "medication_data")
)
```

## Arguments

- object:

  a S4 object of class "metime_analyser".

- which_data:

  character vector of dataset names to convert. If NULL, converts all
  datasets.

- exclude:

  character vector of dataset names to skip. Defaults to
  phenotype/medication data.

## Value

A named list of SummarizedExperiment objects.
