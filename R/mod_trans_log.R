
#' Function to apply log transformation
#' @description Modification(mod) function to log transform dataset/s
#' @examples # example to apply log transformation
#' object <- mod_logtrans(object, which_data="name of the dataset", base=2)
#' @param object An object of class metime_analyser
#' @param which_data Name/s of the dataset to be used
#' @param base base of logarithm to be used
#' @return An object of class metime_analyser with log-transformed dataset/s
#' @seealso [base::log], [mod_trans_zscore]
#' @export
setGeneric("mod_trans_log", function(object, which_data, base) standardGeneric("mod_trans_log"))
setMethod("mod_trans_log", "metime_analyser",function(object, which_data, base=2) {
  #define data to be processed
  data_position <- which(names(object@list_of_data) %in% which_data)
  object@list_of_data[data_position] = lapply(object@list_of_data[data_position] , log, base=base)
  out <- object
  out <- add_function_info(object=out, function_name="mod_trans_log", params=list(which_data=which_data, base=base))
  return(out)
})

