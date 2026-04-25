# Pack all the data into a single object of class "metime_analyser"

Create an object of class "metime_analyser" from a dataset including
data, row_data (id column corresponds to rownames(data)) and col_data
(id colummn corresponds to colnames(data))

## Usage

``` r
get_make_analyser_object(
  data,
  col_data,
  row_data,
  annotations_index = list(),
  name = "set_1",
  results = list()
)
```

## Arguments

- data:

  a data.frame containing data.

- col_data:

  a data.frame containing col_data: id column of col data has to match
  colnames of data.

- row_data:

  a data.frame containing row_data: id column of row data has to match
  rownames of data.

- annotations_index:

  a named list to be filled as list(phenotype="Name or index of the
  file/list", medication="Name or index of the files/list").

- name:

  a character to be assign to the new dataset. Default is set to
  "set_1".

- results:

  a list of existing results. Default set to NULL.

## Value

An object of class metime_analyser
