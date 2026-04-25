# An automated fucntion to calculate GGM from genenet crosssectional version

automated funtion that can be applied on metime_analyser object to
obtain geneNet network along with threshold used

## Usage

``` r
calc_ggm_genenet(
  object,
  which_data,
  threshold,
  all,
  cols_for_meta,
  covariates,
  stratifications,
  name,
  ...
)
```

## Arguments

- object:

  S4 object of cĺass metime_analyser

- which_data:

  a character or a character vector naming the datasets of interest

- threshold:

  type of threshold to be used for extracting significant edges. allowed
  inputs are "li", "FDR", "bonferroni"

- all:

  Logical to extract all edges without any pval correction

- cols_for_meta:

  list of character vector for extracting metadata of metabolites for
  plotting

- covariates:

  covariates to be used for this analysis

- stratifications:

  List to stratify data into a subset. Usage list(name=value)

- name:

  Name of the result

- ...:

  additional arguments for GeneNet

## Value

Network data as a plotter object
