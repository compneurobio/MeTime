# Extending Modules

> One-line summary: add new `mod_*`, `calc_*`, or `meta_*` methods safely.

## Naming conventions

- Use prefixes consistently (`add_`, `mod_`, `calc_`, `meta_`, `get_`, `write_`).

## Skeleton: `mod_*/add_*`


```r

setGeneric("[add|mod]_new_name", function(object, ...) standardGeneric("[add|mod]_new_name"))

setMethod("[add|mod]_new_name", "metime_analyser", function(object, ...) {
            # conditions to check if the arguments if an argument fails return the object as is
            # An example is given here
            if(!all(which_data %in% names(object@list_of_data))) {
                warning("dataset not found in the object. Exiting without making any changes")
                return(object)
            }

            # Add Logic to add/modify elements to the object

            # finally update the functions_applied section of the results
            out <- object
            out <- add_function_info(object=out, function_name="[add|mod]_new_name", params=list(...))
            return(out)
    })

```


## Skeleton: `calc_*`


```r

setGeneric("calc_new_name", function(object, which_data, cols_for_meta=NULL, name="calc_new_name_1", stratifications=NULL, ...) standardGeneric("calc_new_name"))

setMethod("calc_new_name", "metime_analyser", function(object, which_data, cols_for_meta=NULL, name="calc_new_name_1", stratifications=NULL, ...) {
            # conditions to check if the arguments if an argument fails return the object as is
            # An example is given here, of course there might be many more based on arguments.
            if(!all(which_data %in% names(object@list_of_data))) {
                warning("dataset not found in the object. Exiting without making any changes")
                return(object)
            }

            # Updating the name of the result here. Make sure to check the regex based on the name of the function
            if(grep(name, names(object@results)) %>% length() >=1) {
                warning("name of the results was previously used, using a different name")
                index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
                index <- c(0:9)[grep(index, 0:9)+1]
                name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
            }

            # Performing stratification analysis before proceeding with calculation
            data_list <- get_stratified_data(object=object, which_data=which_data, stratifications=stratifications)
            data <- data_list[["data"]]
            row_data <- data_list[["row_data"]]

            # get metadata based on cols_for_meta. These columns will later be used as color or shape aesthetics while plotting the results. Depending on the calculation choose metadata for metabolites or metadata for samples using get_metadata_for_columns() and get_metadata_for_rows() respectively. An example with metabolites is shown here
            metadata <- get_metadata_for_columns(object=object, which_data=which_data, columns=cols_for_meta)

            # Add Logic for calculation on data and row_data generated above here

            # Update results and function information here
            out <- get_make_results(object=object, data = result_here # this is a list, 
                                metadata = metadata, 
                                calc_type = rep("name_of_calc", each=length(result_here)), 
                                calc_info = "More detailed information as a string here",
                                name=name) %>%
                    add_function_info(object=out, function_name="calc_new_name", 
                            params=list(which_data=which_data, ..., cols_for_meta=cols_for_meta, 
                            name=name, stratifications=stratifications))
            return(out)
    })



```

This example is restricted to single dataset. However for allowing multiple datasets at once wrap the code chunk from stratification to add_function_info() part using a loop. Also, make sure that the column names mentioned in stratifications list is consistent across the row_data of the datasets. To see such an example look at calc_conservation_* functions.


## Skeleton: `get_*/write_*`

get_* functions have no strict template as the functions don't manipulate the object itself

## Skeleton: `meta_*`

meta_* functions currently store the results in a different class and currently not extendable
