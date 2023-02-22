#' Function to mutate columns in row_data or col_data
#' @description This function allows you to mutate columns by changing class or ...
#' @param object An S4 object of class metime_analyser
#' @param which_data Dataset of interest. Has to be of length=1
#' @param type either "row" for row_data or "col" for col_data or "data" for data. Set to "data" as default
#' @param ... arguments to pass directly into dplyr::mutate() function.
#' @returns object with mutated col_data and row_data
#' @export
setGeneric("mod_mutate", function(object, which_data, type="data", ...) standardGeneric("mod_mutate"))
setMethod("mod_mutate", "metime_analyser", function(object, which_data, type="data", ...) {
		stopifnot(length(which_data)==1)
		stopifnot(which_data %in% names(object@list_of_data))
		# add_function_info "..." will be in params
		stopifnot(type %in% c("data", "row_data", "col_data"))
		if(type %in% "data") {
			data <- object %>% get_data(which_data=which_data)
			data <- data %>% dplyr::mutate(...)
			object@list_of_data[[which_data]] <- data
		} else if(type %in% "col_data") {
			data <- object %>% get_coldata(which_data=which_data)
			data <- data %>% dplyr::mutate(...)
			object@list_of_col_data[[which_data]] <- data
		} else if(type %in% "row_data") {
			data <- object %>% get_rowdata(which_data=which_data)
			data <- data %>% dplyr::mutate(...)
			object@list_of_row_data[[which_data]] <- data
		}
		return(object)
	})