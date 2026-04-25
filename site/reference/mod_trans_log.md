# Function to apply log transformation

Modification(mod) function to log transform dataset/s

## Usage

``` r
mod_trans_log(object, which_data, base)
```

## Arguments

- object:

  An object of class metime_analyser

- which_data:

  Name/s of the dataset to be used

- base:

  base of logarithm to be used

## Value

An object of class metime_analyser with log-transformed dataset/s

## See also

base::log,
[mod_trans_zscore](https://compneurobio.github.io/MeTime/reference/mod_trans_zscore.md)

## Examples

``` r
# example to apply log transformation
object <- mod_logtrans(object, which_data="name of the dataset", base=2)
#> Error in mod_logtrans(object, which_data = "name of the dataset", base = 2): could not find function "mod_logtrans"
```
