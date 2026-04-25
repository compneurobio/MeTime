# Function to add the clusters obtained from wgcna

Modification(mod) function to get cluster assiginment and generate
eigenmetabolite matrix from a dataset.

## Usage

``` r
mod_trans_eigendata(
  object,
  which_data,
  append,
  results_index = NULL,
  cols_for_meta = NULL,
  baseline = "t0",
  name = "WGCNA_clusters_1",
  ...
)
```

## Arguments

- object:

  An S4 object of class metime_analyser

- which_data:

  character to define which dataset is to be used

- append:

  logical if set to true adds the new data to the object used else
  creates new object

- results_index:

  results_index if clusters were previously calculated else set to
  NULL(default)

- cols_for_meta:

  A list of named character vector to extract col_data to add it for the
  eigendata. will be parsed to get_metadata_for_columns

- baseline:

  a character to define the baseline time point which is used for
  cluster calculation.

- ...:

  arguments for calc_clusters_wgcna. Make sure to set the correct
  baseline value if you are using the function directly

## Value

metime_analyser object with new dataset with eigendata of the
metabolites

## See also

[get_metadata_for_columns](https://compneurobio.github.io/MeTime/reference/get_metadata_for_columns.md),
[calc_clusters_wgcna](https://compneurobio.github.io/MeTime/reference/calc_clusters_wgcna.md)
