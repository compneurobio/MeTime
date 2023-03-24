
#' Function to scale the data
#' @description Modification(mod) Function for scaling datasets 
#' @examples # example to apply scaling
#' object <- mod_trans_zscore(object, which_data=c("dataset1", ...))
#' @param object An object of class metime_analyser
#' @param which_data character vector to define the dataset/s to be used
#' @seealso [base::scale], [mod_trans_log]
#' @return An object of class metime_analyser with scaled which_data
#' @export
setGeneric("mod_trans_zscore", function(object, which_data) standardGeneric("mod_trans_zscore"))
setMethod("mod_trans_zscore", "metime_analyser", function(object, which_data) {
  #sanity checks
  if(!all(which_data %in% names(object@list_of_data))) {
    warning("mod_trans_zscore(): datasets mentioned are not found in the data. Exiting without making any changes.")
    return(object)
  }
  #define data to be processed
  data_position <- which(names(object@list_of_data) %in% which_data)
  data_rownames <- lapply(object@list_of_data[data_position], function(a) rownames(a))
  object@list_of_data[data_position] <- lapply(object@list_of_data[data_position], function(a) {
          a <- apply(a, 2, as.numeric)
          return(a)
    })
  object@list_of_data[data_position] = lapply(object@list_of_data[data_position], scale, center=TRUE, scale=TRUE)
  object@list_of_data[data_position] = lapply(object@list_of_data[data_position], as.data.frame)
  for(i in 1:length(data_position)) {
      rownames(object@list_of_data[[data_position[i]]]) <- data_rownames[[i]] 
  }
  out <- object
  out <- add_function_info(object=out, function_name="mod_trans_zscore", params=list(which_data=which_data))
  return(out)
})


