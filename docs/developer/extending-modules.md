# Extending Modules

> One-line summary: add new `mod_*`, `calc_*`, or `meta_*` methods safely.

## Naming conventions

- Use prefixes consistently (`add_`, `mod_`, `calc_`, `meta_`, `get_`, `write_`).

## Skeleton: `mod_*`

```r
setGeneric("mod_new_name", function(object, ...) standardGeneric("mod_new_name"))
setMethod("mod_new_name", "metime_analyser", function(object, ...) {
  # TODO: validate args
  # TODO: modify object
  # TODO: add function info
  object
})
```

## Skeleton: `calc_*`

```r
setGeneric("calc_new_name", function(object, which_data, ...) standardGeneric("calc_new_name"))
setMethod("calc_new_name", "metime_analyser", function(object, which_data, ...) {
  # TODO: validate inputs
  # TODO: calculate
  # TODO: write results metadata
  object
})
```
