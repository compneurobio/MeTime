#' Function to add information of function added to the data
#' @description Function to add information about the method applied to the dataset
#' @param object S4 object of class metime_analyser
#' @param function_name name of the function used
#' @param params other parameters used wrapped as a list
#' @return object of class metime_analyser with the information of method applied
#' @export 
setGeneric("add_function_info", function(object, function_name, params) standardGeneric("add_function_info"))
setMethod("add_function_info", "metime_analyser", function(object, function_name, params) {
			params <- unlist(params)
			results <- object@results
			if(length(results)==0 | length(results)==1) {
				if(length(results)==0) results[[1]] <- list()
				if(length(results)==1) {
					if(length(grep("calc_|mod_merge_results", names(results[[1]]$functions_applied)))==0) {
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
				if(length(grep("calc_|mod_merge_results", results[[length(results)]]$functions_applied))==0) {
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
					results[[length(results)+1]]$functions_applied[[1]] <- params
					names(results[[length(results) + 1]]$functions_applied)[1] <- function_name
				} 
			}
			object@results <- results
			if(validObject(object)) {
				return(object)
			} else {
				stop("Issue in creating a valid object of class metime_analyser")
			}
	}) 