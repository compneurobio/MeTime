# Function to scale the data

Modification(mod) Function for scaling datasets

## Usage

``` r
mod_trans_zscore(object, which_data)
```

## Arguments

- object:

  An object of class metime_analyser

- which_data:

  character vector to define the dataset/s to be used

## Value

An object of class metime_analyser with scaled which_data

## See also

[base::scale](https://rdrr.io/r/base/scale.html),
[mod_trans_log](https://compneurobio.github.io/MeTime/reference/mod_trans_log.md)

## Examples

``` r
# example to apply scaling
object <- mod_trans_zscore(object, which_data=c("dataset1", ...))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'mod_trans_zscore': object 'object' not found
```
