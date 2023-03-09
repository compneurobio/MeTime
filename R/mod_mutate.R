#' Function to mutate columns in row_data or col_data
#' @description This function allows you to mutate columns by changing class or ...
#' @param object An S4 object of class metime_analyser
#' @param which_data Dataset of interest. Has to be of length=1
#' @param type either "row" for row_data or "col" for col_data or "data" for data and "results" for results. 
#' Set to "data" as default
#' @param results_index Index to define the results of interest. Can be character and numeric
#' @param ... arguments to pass directly into dplyr::mutate() function.
#' @returns object with mutated col_data and row_data
#' @export
setGeneric("mod_mutate", function(object, which_data, type="data", results_index=NULL, ...) standardGeneric("mod_mutate"))
setMethod("mod_mutate", "metime_analyser", function(object, which_data, type="data", results_index=NULL, ...) {
		stopifnot(length(which_data)==1)
		stopifnot(which_data %in% names(object@list_of_data))
		stopifnot(type %in% c("data", "row_data", "col_data", "results"))
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
			if(is.null(results_index)) {
				warning("Results index is not specified, exiting without making any changes")
				return(object)
			}
			if(grep("ggm|network", object@results$information$calc_type) %>% length() >=1) {
				warning("mutations and filters are not possible for networks. Exiting without making any changes")
				return(object)
			}
			results <- object@results[[results_index]]
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