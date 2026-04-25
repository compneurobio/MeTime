# Add phenotypic measurements that are not added to the row_data of the dataset

A method applied on the s4 object of class "metime_analyser" to add all
those datapoints that are present in phenotype dataframe but not in
row_data comes with the feature of updating those data points measured
only at screening to all datapoints and then adding it to row_data

## Usage

``` r
add_distribution_vars_to_rows(
  object,
  screening_vars,
  distribution_vars,
  which_data
)
```

## Arguments

- object:

  An object of class metime_analyser

- screening_vars:

  A character vector to define the vars that are to be updated as per
  add_screening_vars() else set it to NULL(Default).

- distribution_vars:

  A character naming the vars of interest

- which_data:

  dataset to which the information is to be added(only 1 can be used at
  a time)

## Value

object of class metime_analyser with phenotype data added to row data

## See also

[add_screening_vars](https://compneurobio.github.io/MeTime/reference/add_screening_vars.md)

## Examples

``` r
object <- add_distribution_vars_to_rows(object=data, screening_vars=c("var1", "var2"), 
    distribution_vars=c("var1", "var2", "var3"), which_data="dataset1")
#> Error in (function (classes, fdef, mtable) {    methods <- .findInheritedMethods(classes, fdef, mtable)    if (length(methods) == 1L)         return(methods[[1L]])    else if (length(methods) == 0L) {        cnames <- paste0("\"", vapply(classes, as.character,             ""), "\"", collapse = ", ")        stop(gettextf("unable to find an inherited method for function %s for signature %s",             sQuote(fdef@generic), sQuote(cnames)), domain = NA)    }    else stop("Internal error in finding inherited methods; didn't return a unique method",         domain = NA)})(list("function"), new("standardGeneric", .Data = function (object,     screening_vars, distribution_vars, which_data) standardGeneric("add_distribution_vars_to_rows"), generic = structure("add_distribution_vars_to_rows", package = "MeTime"),     package = "MeTime", group = list(), valueClass = character(0),     signature = c("object", "screening_vars", "distribution_vars",     "which_data"), default = NULL, skeleton = (function (object,         screening_vars, distribution_vars, which_data)     stop(gettextf("invalid call in method dispatch to '%s' (no default method)",         "add_distribution_vars_to_rows"), domain = NA))(object,         screening_vars, distribution_vars, which_data)), <environment>): unable to find an inherited method for function ‘add_distribution_vars_to_rows’ for signature ‘"function"’
```
