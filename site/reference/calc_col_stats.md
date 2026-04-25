# Check normality of metabolites and append results to col_data

A method applied to data within the S4 object of class "metime_analyser"
to check normality of the metabolites and generate a results section for
easy visualization,

## Usage

``` r
calc_col_stats(
  object,
  which_data,
  method = c("shapiro", "ks"),
  cols_for_meta,
  name = "calc_col_stats_1"
)
```

## Arguments

- object:

  a S4 object of the class "metime_analyzer".

- which_data:

  a character to define which dataset is to be used.

- method:

  a character vector specifying the method to check for normality. By
  default all available options are chosen: c("shapiro", "ks")". See
  details for more information.

- cols_for_meta:

  a list of named character vectors for obtaining metainformation of
  metabolites

- name:

  name of the results. Default is set to "calc_col_stats_1"

## Value

a S4 object of class "metime_analyser" with results appended to the
col_data of which_data

## Details

Shapiro: Performs the Shapiro-Wilk test of normality (normal
distribution if shapiro_pval \> 0.05). This test is based on correlation
between the data and the corresponding normal scores. Ks: Performs a
one-sample Kolmogorov-Smirnov test (normal distribution if ks_pval \>
0.05). The ks.test compares the metabolite values against pnorm.

## See also

[shapiro.test](https://rdrr.io/r/stats/shapiro.test.html)
[ks.test](https://rdrr.io/r/stats/ks.test.html)
[get_metadata_for_columns](https://compneurobio.github.io/MeTime/reference/get_metadata_for_columns.md)
