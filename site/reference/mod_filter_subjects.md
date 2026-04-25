# Filter subjects from datasets

Removes samples that are not in the specified subjects.

## Usage

``` r
mod_filter_subjects(object, subjects, which_data = NULL)
```

## Arguments

- object:

  a S4 object of class "metime_analyser".

- subjects:

  character vector of subjects to keep.

- which_data:

  character vector of dataset names to filter. If NULL, filters all
  datasets.

## Value

A modified metime_analyser object.
