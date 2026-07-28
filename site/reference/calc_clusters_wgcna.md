# Calculate cluster assignment from WGCNA

calculate cluster assignment to specified col_data from Weighted
Correlation Network Analysis

## Usage

``` r
calc_clusters_wgcna(
  object,
  which_data,
  baseline = "t0",
  cols_for_meta = NULL,
  name = "WGCNA_clusters_1",
  soft_power = NULL,
  ...
)
```

## Arguments

- object:

  a S4 object of the class metime_analyzer.

- which_data:

  a character to define which dataset is to be used.

- baseline:

  a character to define the baseline time point which is used for
  cluster calculation.

- cols_for_meta:

  A list of named character vector of length equal to which_data.
  Example: list(dataset1=c(new_name1="colname", new_name2="colname"))
  Default is set to NULL

- name:

  Name of the results. Default is set to "WGCNA_clusters_1"

- soft_power:

  Numeric soft-thresholding power to use. When `NULL`, the first
  candidate power with a truncated R-squared greater than 0.8 is used.
  If no candidate meets that threshold, the function warns and returns
  without changing `object`. Supply a numeric value to override the
  automatic choice.

- ...:

  multiple parameters separated by commas passed to cutreeDynamic (R
  package: dynamicTreeCut). Parameters include: minClusterSize,
  pamDendroRespect ...

## Value

a S4 object of class "metime_analyser" with cluster information appended
to col_data of which_data

## Details

Based on the method described in the WGCNA tutorials for step-by-step
network construction and module detection. The idea is to find modules
of metabolites at baseline. This method can be found in detail here
[WGCNA: an R package for weighted correlation network
analysis](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)

## See also

[dynamicTreeCut::cutreeDynamic](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamic.html),
[get_metadata_for_columns](https://compneurobio.github.io/MeTime/reference/get_metadata_for_columns.md),
[mod_trans_eigendata](https://compneurobio.github.io/MeTime/reference/mod_trans_eigendata.md)
