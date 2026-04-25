# Function to extract analyser object data into a csv

extracts information from analyser object and saves it as a csv

## Usage

``` r
write_data(object, which_data, type)
```

## Arguments

- object:

  An object of class metime_plotter

- which_data:

  Character to specify the dataset

- type:

  which type of output file. Can be "csv", "tsv" and "xlsx"

## Value

saves the data in the working directory as a csv and returns nothing

## See also

[write_results](https://compneurobio.github.io/MeTime/reference/write_results.md),
[write_report](https://compneurobio.github.io/MeTime/reference/write_report.md)

## Examples
