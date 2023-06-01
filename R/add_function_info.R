#' Add information of function added to the data
#' @description Add information on the method applied to the metime_analyzer object.
#' @param object a S4 object of the class "metime_analyzer".
#' @param function_name a character of the name of the function that was used.
#' @param params list of parameters that were used in function.
#' @return object of class metime_analyser with the information of method applied
#' @export 
setGeneric("add_function_info", function(object, function_name, params) standardGeneric("add_function_info"))
setMethod("add_function_info", "metime_analyser", function(object, function_name, params) {
			my_list <- purrr::map(params, ~rlang::enquos(.x)) %>% unlist(recursive=FALSE)
			out <- lapply(seq_along(my_list),function(x) paste0(names(my_list)[x],as.character(my_list[x]))) %>% 
                  		unlist() %>% 
    						gsub(pattern="~", replacement=" = ") %>% 
    						paste0(collapse=",") %>% 
    						paste0(function_name, "(",.,")") %>%
    						noquote()
			results <- object@results
			if(length(results)==0 | length(results)==1) {
				if(length(results)==0) results[[1]] <- list()
				if(length(results)==1) {
					if(length(grep("calc_|mod_merge_results|add_result|meta_", results[[1]]$functions_applied))==0) {
						if(length(results[[1]]$functions_applied) >=1) {
							results[[1]]$functions_applied[length(results[[1]]$functions_applied)+1] <- out
						} else {
							results[[1]]$functions_applied[1] <- out
						}
					} else {
						results[[2]] <- list()
						results[[2]]$functions_applied[1] <- out
					}
				}
			} else {
				if(length(grep("calc_|mod_merge_results|add_result|meta_", results[[length(results)]]$functions_applied))==0) {
					if(length(results[[length(results)]]$functions_applied) >=1) {
						results[[length(results)]]$functions_applied[length(results[[length(results)]]$functions_applied)+1] <- out
					} else {
						results[[length(results)]]$functions_applied[1] <- out
					}
				} else {
					results[[length(results) + 1]] <- list()
					results[[length(results)]]$functions_applied[1] <- out
				} 
			}
			object@results <- results
			if(validObject(object)) {
				return(object)
			} else {
				stop("Issue in creating a valid object of class metime_analyser")
			}
	}) 