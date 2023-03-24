#' Function to rename columns in the datasets and results of analyser object
#' @description Function to rename columns of data, row_data, col_data and plot_data of the analyser object.
#' Arguments will be passed directly into dplyr::rename() function
#' @param object An S4 object of class metime_analyser
#' @param which_data character or numeric input of length 1 to define the index of dataset or results
#' @param type character input to define the type of data to be manipulated. Accepted inputs are "row_data", "col_data", "data"
#' and "results". However renamed results will be returned to the user as a list of results and will not return the full analyser
#' object
#' @param ... arguments to pass to dplyr::rename to change the name of the functions. Example new_colname = "old_colname"
#' @seealso [mod_filter] [mod_mutate]
#' @export

setGeneric("mod_rename", function(object, which_data, type="data", ...) standardGeneric("mod_rename"))
setMethod("mod_rename", "metime_analyser", function(object, which_data, type="data", ...) {
		if(length(which_data)!=1) {
			warning("length of which_data is not 1, exiting without making any changes")
			return(object)
		}
		if(!type %in% c("data", "row_data", "col_data", "results")) {
			warning("type is unknown. Please check the accepted values and exiting without making any changes")
			return(object)
		}
		if(type %in% "data") {
			object@list_of_data[[which_data]] <- object@list_of_data[[which_data]] %>% dplyr::rename(...)
			new_names <- list(...)
			old_names <- unname(new_names) %>% unlist()
			new_names <- names(new_names) %>% as.character()
			object@list_of_col_data[[which_data]]$id <- ifelse(object@list_of_col_data[[which_data]]$id %>% old_names,
							new_names,
							object@list_of_col_data[[which_data]]$id)
			rownames(object@list_of_col_data[[which_data]]) <- object@list_of_col_data[[which_data]]$id
		} else if(type %in% "row_data") {
			object@list_of_row_data[[which_data]] <- object@list_of_row_data[[which_data]] %>% dplyr::rename(...)
		} else if(type %in% "col_data") {
			object@list_of_col_data[[which_data]] <- object@list_of_col_data[[which_data]] %>% dplyr::rename(...)
		} else if(type %in% results) {
			plot_data <- lapply(seq_along(object@results[[which_data]]$plot_data), function(x) {
						object@results[[which_data]]$plot_data[[x]] %>% dplyr::rename(...) %>% return()
				})
			return(plot_data)
		}
		exprs <- as.list(substitute(list(...))[-1])
		str <- paste0(names(exprs), "=", sapply(exprs, function(x) {
		   if (is.character(x)) paste0('"', x, '"')
		   else as.character(x)
		 }), collapse = ", ")
		object <- object %>% add_function_info(function_name="mod_rename", params=list(which_data=which_data, type=type, mutations=str))
		return(object)

	})