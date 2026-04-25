# Function to stratify before calculation

get function to perform stratification of a dataset in calculations.
Used only in calc\_ functions.

## Usage

``` r
get_stratified_data(object, which_data, stratifications)
```

## Arguments

- object:

  An S4 object of class metime_analyser

- which_data:

  character vector to define Dataset/s to be used

- stratifications:

  list of variables and their values to stratified - list(name=value).
  example of stratification with time
  stratifications=list(time=c("timepoint1", "timepoint2", ...)) then
  only these timepoints will be considered in the analysis

## Value

data_list that contains col_data, row_data and data with stratifications
aforementioned. Access them as shown: data=data_list\["data"\];
col_data=data_list\["col_data"\]; row_data=data_list\["row_data"\]
