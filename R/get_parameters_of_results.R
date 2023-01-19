#' Get parameters of the results
#' @description Function to get parameters and functions applied to obtain results
#' @param object An S4 object of class metime_analyser
#' @param results_index name or index to get to the results of interest
#' @return a dataframes with functions(rownames) and parameters(colnames)
#' @export
setGeneric("get_parameters_of_results", function(object, results_index) standardGeneric("get_parameters_of_results"))
setMethod("get_parameters_of_results", "metime_analyser", function(object, results_index) {
		if(class(results_index) %in% "character") stopifnot(results_index %in% names(object@results))
		if(class(results_index) %in% c("numeric", "integer")) stopifnot(results_index < length(object@results))
		parameters <- object@results[[results_index]][["functions"]]
		functions <- names(parameters)
		out <- lapply(seq_along(parameters), function(i) {
				func <- parameters[[i]]
				func <- enframe(func) %>% unnest() %>% group_by(name) %>% mutate(value=paste(value, collapse=",")) %>%
						distinct() %>% t() %>% as.data.frame()
				colnames(func) <- func[1,]
				func <- func[-1, ]
				return(func)
			}) %>% do.call(what=plyr::rbind.fill)
		rownames(out) <- functions
		return(out)
	})