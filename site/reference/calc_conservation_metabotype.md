# Function to calculate metabotype conservation index

Method applied on the object metime_analyser to calculate the metabotype
conservation index

## Usage

``` r
calc_conservation_metabotype(
  object,
  which_data,
  verbose = F,
  cols_for_meta,
  stratifications,
  name = "calc_conservation_metabotype_1"
)
```

## Arguments

- object:

  An object of class metime_analyser

- which_data:

  Name of the dataset to be used

- verbose:

  Information provided on steps being processed

- cols_for_meta:

  Character vector to define column names that are to be used for
  plotting purposes

- stratifications:

  List to stratify data into a subset. Usage list(name=value)

- name:

  character vector to define the results. Should be equal to length of
  which_data

## Value

List of conservation index results

## Examples

``` r
#calculating metabotype_conservation_index 
out <- calc_metabotype_conservation(object=metime_analyser_object, which_data="Name of the dataset")
#> Error in calc_metabotype_conservation(object = metime_analyser_object,     which_data = "Name of the dataset"): could not find function "calc_metabotype_conservation"
```
