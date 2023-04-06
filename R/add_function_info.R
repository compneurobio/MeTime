#' Add information of function added to the data
#' @description Add information on the method applied to the metime_analyzer object.
#' @param object a S4 object of the class "metime_analyzer".
#' @param function_name a character of the name of the function that was used.
#' @param params list of parameters that were used in function.
#' @return object of class metime_analyser with the information of method applied
#' @export 
setGeneric("add_function_info", function(object, function_name, params) standardGeneric("add_function_info"))
setMethod("add_function_info", "metime_analyser", function(object, function_name, params) {
			params <- unlist(params)
			results <- object@results
			if(length(results)==0 | length(results)==1) {
				if(length(results)==0) results[[1]] <- list()
				if(length(results)==1) {
					if(length(grep("calc_|mod_merge_results|add_result|meta_", names(results[[1]]$functions_applied)))==0) {
						if(length(results[[1]]$functions_applied) >=1) {
							names <- names(results[[1]]$functions_applied)
							results[[1]]$functions_applied[[length(results[[1]]$functions_applied)+1]] <- params
							names(results[[1]]$functions_applied) <- c(names, function_name)
						} else {
							results[[1]]$functions_applied[[1]] <- params
							names(results[[1]]$functions_applied)[1] <- function_name
						}
					} else {
						results[[2]] <- list()
						results[[2]]$functions_applied[[1]] <- params
						names(results[[2]]$functions_applied)[1] <- function_name
					}
				}
			} else {
				if(length(grep("calc_|mod_merge_results|add_result|meta_", names(results[[length(results)]]$functions_applied)))==0) {
					if(length(results[[length(results)]]$functions_applied) >=1) {
						names <- names(results[[length(results)]]$functions_applied)
						results[[length(results)]]$functions_applied[[length(results[[length(results)]]$functions_applied)+1]] <- params
						names(results[[length(results)]]$functions_applied) <- c(names, function_name)
					} else {
						results[[length(results)]]$functions_applied[[1]] <- params
						names(results[[length(results)]]$functions_applied)[1] <- function_name 
					}
				} else {
					results[[length(results) + 1]] <- list()
					results[[length(results)]]$functions_applied[[1]] <- params
					names(results[[length(results)]]$functions_applied)[1] <- function_name
				} 
			}
			object@results <- results
			if(validObject(object)) {
				return(object)
			} else {
				stop("Issue in creating a valid object of class metime_analyser")
			}
	}) 