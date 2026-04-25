# Function to calculate dimensionality reduction methods such as tsne, umap and pca.

A method to apply on s4 object of class metime_analyse in order to
obtain information after dimensionality reduction on a dataset/s

## Usage

``` r
calc_dimensionality_reduction_metabs(
  object,
  which_data,
  type,
  cols_for_meta,
  stratifications = list(),
  name,
  ...
)
```

## Arguments

- object:

  An object of class metime_analyser

- which_data:

  a character vector - Names of the dataset from which the samples will
  be extracted

- type:

  type of the dimensionality reduction method to be applied. Accepted
  inputs are "UMAP", "tSNE", "PCA"

- cols_for_meta:

  A list of Character vectors to define columns names that are to be
  used for plotting purposes

- stratifications:

  List to stratify data into a subset. Usage list(name=value)

- name:

  name of the results that you want to have

- ...:

  additional arguments that can be passed on to prcomp(), M3C::tsne()
  and umap::umap()

## Value

a list with two plotter objects containing the dimensionality reduction
information that can be parsed into plotting function

## Examples

``` r
#calculate PCA
pca <- calc_dimensionality_reduction(object=metime_analyser_object, which_data="name/s of the dataset/s", type="PCA")
#> Error in calc_dimensionality_reduction(object = metime_analyser_object,     which_data = "name/s of the dataset/s", type = "PCA"): could not find function "calc_dimensionality_reduction"
#calculate UMAP
pca <- calc_dimensionality_reduction(object=metime_analyser_object, which_data="name/s of the dataset/s", type="UMAP")
#> Error in calc_dimensionality_reduction(object = metime_analyser_object,     which_data = "name/s of the dataset/s", type = "UMAP"): could not find function "calc_dimensionality_reduction"
#calculate tSNE
pca <- calc_dimensionality_reduction(object=metime_analyser_object, which_data="name/s of the dataset/s", type="tSNE")
#> Error in calc_dimensionality_reduction(object = metime_analyser_object,     which_data = "name/s of the dataset/s", type = "tSNE"): could not find function "calc_dimensionality_reduction"
```
