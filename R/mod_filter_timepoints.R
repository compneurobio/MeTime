
#' Functions for selecting time points
#' @description a method applied onto class metime_analyser in order to extract timepoints of interest from a dataset
#' @examples #example to use this function
#' object <- mod_filter_tp(object, timepoints=c("t0","t12","t24"), full=TRUE, which_data="Name of the dataset")
#' @param object An object of class metime_analyser
#' @param timepoints time points to be selected. 
#' @param which_data Name of the dataset to be used
#' @param complete if TRUE subjects are only selected if measured in all selected time points
#' @return An object of class metime_analyser with processed data
#' @export 

setGeneric("mod_filter_timepoints", function(object, timepoints, complete, which_data) standardGeneric("mod_filter_timepoints"))
setMethod("mod_filter_timepoints", "metime_analyser", function(object, timepoints, complete=TRUE, which_data) {
  # define data to be processed
  data_position <- which(names(object@list_of_data) %in% which_data)
  
  test <- lapply(data_position, function(i) {
    keep_id <- object@list_of_row_data[[i]] %>% 
      dplyr::select(id, time, subject) %>%  
      dplyr::filter(time %in% timepoints)
    if(full) {
      full_rid <- keep_id %>% 
        dplyr::count(subject) %>% 
        dplyr::filter(n==length(timepoints))
      keep_id <- keep_id %>% 
        dplyr::filter(subject %in% full_rid$subject)
      object@list_of_row_data[[i]] = object@list_of_row_data[[i]] %>% 
        dplyr::filter(id %in% keep_id$id)
      object@list_of_data[[i]] = object@list_of_data[[i]][keep_id$id,]
    } else {
      object@list_of_row_data[[i]] = object@list_of_row_data[[i]] %>% 
        dplyr::filter(id %in% keep_id$id)
      object@list_of_data[[i]] = object@list_of_data[[i]][keep_id$id,]
    }
    return(NULL)
  })
  out <- object
  out <- add_function_info(object=out, function_name="mod_filter_timepoints",
      params=list(timepoints=timepoints, complete=complete, which_data=which_data))
  return(out)
})


