#' Function to add information of function added to the data
#' @decription Function to add information about the method applied to the dataset
#' @param object S4 object of class metime_analyser
#' @param results_index dataset/s to which function was applied
#' @param function_name name of the function used
#' @param params other parameters used wrapped as a list
#' @return object of class metime_analyser with the information of method applied
#' @export 
setGeneric("add_function_info", function(object, function_name, params) standardGeneric("add_function_info"))
setMethod("add_function_info", "metime_analyser", function(object, function_name, params) {
			params <- unlist(params)
			results <- object@results
			if(length(results)==0 | length(results)==1) {
				if(length(grep("calc_", names(results[[1]][["functions"]])))==0) {
					if(length(results[[1]][["functions"]]) >=1) {
						names <- names(results[[1]][["functions"]])
						results[[1]][["functions"]][[length(results[[1]][["functions"]])+1]] <- params
						names(results[[1]][["function"]]) <- c(names, function_name)
					} else {
						results[[1]][["functions"]][[1]] <- params
						names(results[[1]][["functions"]]) <- function_name
					}
				} else {
					results[[2]] <- list()
					results[[2]][["functions"]][[1]] <- params
					names(results[[2]][["functions"]]) <- function_name
				}
			} else {
				if(length(grep("calc_", results[[length(results)]][["functions"]]))==0) {
					if(length(results[[length(results)]][["function"]]) >=1) {
						names <- names(results[[length(results)]][["functions"]])
						results[[length(results)]][["functions"]][[length(results[[length(results)]][["functions"]])+1]] <- params
						names(results[[length(results)]][["functions"]]) <- c(names, function_name)
					} else {
						results[[length(results)]][["functions"]][[1]] <- params
					}
				} else {
					results[[length(results) + 1]] <- list()
					results[[length(results)+1]][["functions"]][[1]] <- params
					names(results[[length(results) + 1]]) <- function_name
				} 
			}
			object@results <- results
			if(validObject(object)) {
				return(object)
			} else {
				stop("Issue in creating a valid object of class metime_analyser")
			}
	}) 