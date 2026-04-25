# Function to calculate dimensionality reduction methods such as tsne, umap and pca.

A method to apply on s4 object of class metime_analyse in order to
obtain information after dimensionality reduction on a dataset/s

## Usage

``` r
calc_dimensionality_reduction_samples(
  object,
  which_data,
  type,
  cols_for_meta,
  stratifications = NULL,
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

  A Character vector to define columns names that are to be used for
  plotting purposes

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
pca <- calc_dimensionality_reduction_samples(object=metime_analyser_object, 
 which_data="name/s of the dataset/s", type="PCA", cols_for_meta=c(), stratifications=NULL, name="name of the results")
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'calc_dimensionality_reduction_samples': object 'metime_analyser_object' not found
```
