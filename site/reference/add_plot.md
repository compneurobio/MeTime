# Generates a report from an metime_analyser object

Write an HTML or PDF report that summarizes one or more
"me_time_analyser" objects and display all results

## Usage

``` r
add_plot(
  object,
  result_item = NULL,
  plot_data = NULL,
  calc_type = NA,
  calc_info = NA,
  result_index = "result_index",
  result_type = "result_type",
  result_subtype = "result_subtype",
  result_title = "result_title"
)
```

## Arguments

- object:

  a S4 object of class metime_analyser

- result_item:

  a plot or data table to be added as a result.

- plot_data:

  a dataframe corresponding to the result item.

- calc_type:

  a character defining the calc type of an object that has been added.
  Default set to NA.

- calc_info:

  a character defining the calc info of an object that has been added.
  Default set to NA.

- result_index:

  a character defining the result index of an object that has been
  added. Default set to 'result_index'

- result_type:

  a character defining the result type of an object that has been added.
  Default set to 'result_type'

- result_subtype:

  a character defining the result subtype of an object that has been
  added. Default set to 'result_subtype'

- result_title:

  a character defining the result title of an object that has been
  added. Default set to 'result_title'

## Value

a S4 object of class metime_analyser
