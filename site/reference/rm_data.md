# Function to remove datasets or results from the metime analyser object

An S4 method to remove unwanted data

## Usage

``` r
rm_data(object, which_data, type = "dataset")
```

## Arguments

- object:

  An S4 object of class metime_analyser

- which_data:

  index or name of the dataset to be removed or result

- type:

  "result" or "dataset" based on the type to be removed

## Value

object with results/dataset removed
