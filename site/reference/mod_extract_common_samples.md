# Function to extract common samples across all datasets and store them only

Modification(mod) function applied on object of class metime_analyser to
extract common samples across datasets.

## Usage

``` r
mod_extract_common_samples(object)
```

## Arguments

- object:

  An object of class metime_anaylser

## Value

metime_analyser object with only common samples across all datasets
present in the object parsed

## Examples

``` r
# extracting common samples across all datasets
new_object_with_only_common_samples <- mod_extract_common_samples(object=metime_analyser_object)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'mod_extract_common_samples': object 'metime_analyser_object' not found
```
