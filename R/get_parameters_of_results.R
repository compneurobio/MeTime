#' Get parameters of the results
#' @description Get parameters and functions applied to obtain results.
#' @param object An S4 object of class metime_analyser.
#' @param results_index a character or numeric value to be used as the result index.
#' @return a list with all the function names and the values of the arguments passed into it
#' @export
setGeneric("get_parameters_of_results", function(object, results_index) standardGeneric("get_parameters_of_results"))
setMethod("get_parameters_of_results", "metime_analyser", function(object, results_index) {
		if(class(results_index) %in% "character") stopifnot(results_index %in% names(object@results))
		if(class(results_index) %in% c("numeric", "integer")) stopifnot(results_index < length(object@results))
		my_strings <- object@results[[results_index]][["functions_applied"]]
		extract_arguments <- function(str) {
  			arg_string <- stringr::str_match(str, "[(](.*)[)]")[, 2]
  			eval(parse(text = paste0("list(", arg_string, ")")))
		}
		extract_function_names <- function(str) {
 			match <- str_match(str, '^(\\w+)')
  			if (!is.na(match[1,2])) {
    			return(match[1,2])
  			}
		}
		function_names <- lapply(my_strings, extract_function_names)
		out <- lapply(my_strings, extract_arguments)
		names(out) <- function_names
		return(out)
	})