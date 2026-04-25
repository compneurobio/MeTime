# Generates a report from an metime_analyser object

Write an HTML or PDF report that summarizes one or more
"me_time_analyser" objects and display all results

## Usage

``` r
write_report(
  object,
  title = NULL,
  file = NULL,
  write_results = F,
  device = "html",
  interactive = F,
  author = NULL
)
```

## Arguments

- object:

  a S4 object of class metime_analyser

- title:

  a character string used as title of the report. By default is set to
  metime + system time.

- file:

  a character string used as file name of the output report and table.
  By default is set to metime + system time.

- write_results:

  a logical indicating whether you want to save the plot data as an xlsx
  file. If TRUE an xlsx file with the file name specified by the title
  will be saved.

- device:

  a character string specifying the format of the report. By default is
  set to html. Other options include pdf.

- interactive:

  a logical indicating whether plots are interactive - only possible for
  html. Default set to FALSE

- author:

  a character specifying the name of the author

## Value

Saves the report as html/pdf

## See also

write_report,
[write_results](https://compneurobio.github.io/MeTime/reference/write_results.md)
