#' Get parameters of the results
#' @description Get parameters and functions applied to obtain results.
#' @param object An S4 object of class metime_analyser.
#' @param index a character or numeric value to be used as the result index.
#' @return a dataframe with functions(rownames) and parameters(colnames).
#' @export
setGeneric("get_parameters_of_results", function(object, index) standardGeneric("get_parameters_of_results"))
setMethod("get_parameters_of_results", "metime_analyser", function(object, index) {
		if(class(index) %in% "character") stopifnot(index %in% names(object@results))
		if(class(index) %in% c("numeric", "integer")) stopifnot(index < length(object@results))
		parameters <- object@results[[index]][["functions_applied"]]
		functions <- names(parameters)
		out <- lapply(seq_along(parameters), function(i) {
				func <- parameters[[i]]
				func <- enframe(func) %>% unnest() %>% group_by(name) %>% mutate(value=paste(value, collapse=",")) %>%
						distinct() %>% t() %>% as.data.frame()
				colnames(func) <- func[1,]
				func <- func[-1, ]
				return(as.data.frame(func))
			}) %>% do.call(what=plyr::rbind.fill)
		rownames(out) <- functions
		return(out)
	})