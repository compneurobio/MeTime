#' Function to get results of a particular analysis
#' @description function to get results from S4 object
#' @param object An S4 object of class metime_analyser
#' @param which_result character/index to define the result of interest
#' @returns list with results data list and information regarding 
#' the functions applied and type of calculation
#' @export
setGeneric("get_results", function(object, which_result) standardGeneric("get_results"))
setMethod("get_results", "metime_analyser", function(object, which_result) {
			stopifnot(is.null(which_result))
			return(object@results[[which_result]])
	})