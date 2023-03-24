#' Get results element by index
#' @description Get results from S4 object of class metime_analyser by a specified index. Print object in console to see all available indices.
#' @param object a S4 object of the class "metime_analyzer".
#' @param index a character or numeric value defining which result should be 
#' @returns A list with result elements that 
#' the functions applied and type of calculation
#' @export
setGeneric("get_results", function(object, index) standardGeneric("get_results"))
setMethod("get_results", "metime_analyser", function(object, index) {
  out <- object@results[[index]]
			stopifnot(is.null(which_result))
			return(object@results[[which_result]])
	})