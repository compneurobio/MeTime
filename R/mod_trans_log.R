
#' Function to apply log transformation
#' @description Function to log transform data
#' @examples # example to apply log transformation
#' object <- mod_logtrans(object, which_data="name of the dataset", base=2)
#' @param object An object of class metime_analyser
#' @param which_data Name of the dataset to be used
#' @param base base of log to be used
#' @return An object of class metime_analyser with processed data
#' @export
setGeneric("mod_trans_log", function(object, which_data, base) standardGeneric("mod_trans_log"))
setMethod("mod_trans_log", "metime_analyser",function(object, which_data, base=2) {
  #define data to be processed
  data_position <- which(names(object@list_of_data) %in% which_data)
  object@list_of_data[data_position] = lapply(object@list_of_data[data_position] , log, base=base)
  out <- object
  return(out)
})

