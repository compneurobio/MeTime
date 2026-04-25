# Function to Convert S4 object of class metime_analyser to an S3 object(list) with same architecture

Converter function to be applied onto metime_analyse object to convert
into a standard list of S3 type. For example list_of_data can now be
accessed with '\$' instead of '@'

## Usage

``` r
mod_convert_s4_to_s3(object)
```

## Arguments

- object:

  An object of class metime_analyser

## Value

A list with the same data and structure as a metime_analyser

## Examples

``` r
# convert S4 object to a list
s3_list <- mod_convert_s4_to_s3(object=metime_analyser_object)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'mod_convert_s4_to_s3': object 'metime_analyser_object' not found
```
