#' Function to mutate columns in row_data or col_data
#' @description This function allows you to mutate columns by changing class or ...
#' @param object An S4 object of class metime_analyser
#' @param which_data character or numeric input to define Dataset/result of interest. Has to be of length=1
#' @param type either "row" for row_data or "col" for col_data or "data" for data and "results" for results. 
#' Set to "data" as default
#' @param ... arguments to pass directly into dplyr::mutate() function.
#' @returns object with mutated data
#' @export
setGeneric("mod_mutate", function(object, which_data, type="data", ...) standardGeneric("mod_mutate"))
setMethod("mod_mutate", "metime_analyser", function(object, which_data, type="data", ...) {
		if(length(which_data)!=1) {
			warning("length of which_data is not 1, exiting without making any changes")
			return(object)
		}
		if(!type %in% c("data", "row_data", "col_data", "results")) {
			warning("type is unknown. Please check the accepted values and exiting without making any changes")
			return(object)
		}
		if(type %in% "data") {
			data <- object %>% get_data(which_data=which_data)
			data <- data %>% dplyr::mutate(...)
			object@list_of_data[[which_data]] <- data
			#adjust row data accordingly
    		object@list_of_row_data[[which_data]] <- object %>% 
      			get_rowdata(which_data=which_data) %>%  
      			dplyr::filter(id %in% rownames(object@list_of_data[[which_data]]))
		} else if(type %in% "col_data") {
			data <- object %>% get_coldata(which_data=which_data)
			data <- data %>% dplyr::mutate(...)
			object@list_of_col_data[[which_data]] <- data
			 #adjust data accordingly
   			 object@list_of_data[[which_data]] <- object %>% 
      				get_data(which_data=which_data) %>% 
      				dplyr::select(any_of(object@list_of_col_data[[which_data]]$id))
		} else if(type %in% "row_data") {
			data <- object %>% get_rowdata(which_data=which_data)
			data <- data %>% dplyr::mutate(...)
			object@list_of_row_data[[which_data]] <- data
			#adjust data accordingly
    		object@list_of_data[[which_data]] <- object %>% 
      			get_data(which_data=which_data) %>%
      			dplyr::mutate(id = rownames(.[])) %>% 
      			dplyr::filter(id %in% object@list_of_row_data[[which_data]]$id) %>% 
      			dplyr::select(-id)
		} else if(type %in% "results") {
			#if(is.null(results_index)) {
			#	warning("Results index is not specified, exiting without making any changes")
			#	return(object)
			#}
			if(grep("ggm|network", object@results[[which_data]]$information$calc_type) %>% length() >=1) {
				warning("mutations and filters are not possible for networks. Exiting without making any changes")
				return(object)
			}
			results <- object@results[[which_data]]
			results$plot_data <- lapply(seq_along(results$plot_data), function(x) {
						results$plot_data[[x]] <- results$plot_data[[x]] %>% dplyr::mutate(...)
						return(results$plot_data[[x]])
				})
			return(results$plot_data)
		}
		exprs <- as.list(substitute(list(...))[-1])
		str <- paste0(names(exprs), "=", sapply(exprs, function(x) {
		   if (is.character(x)) paste0('"', x, '"')
		   else as.character(x)
		 }), collapse = ", ")
		object <- object %>% add_function_info(function_name="mod_mutate", params=list(which_data=which_data, type=type, mutations=str))
		return(object)
	})