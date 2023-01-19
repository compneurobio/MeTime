
#' Function to remove duplicates
#' @description Function to remove duplicates from the analyser object
#' @param object An S4 object of class metime_analyser
#' @return object after removing duplicated data
#' @export
setGeneric("mod_remove_duplicates", function(object) standardGeneric("mod_remove_duplicates"))
setMethod("mod_remove_duplicates", "metime_analyser", function(object) {
		out <- lapply(names(object@list_of_data), function(x) {
			data <- object@list_of_data[[x]]
			rowdata <- object@list_of_row_data[[x]]
			coldata <- object@list_of_col_data[[x]]
			data <- data[!duplicated(data), ]
			rowdata <- rowdata[!duplicated(rowdata), ]
			coldata <- coldata[!duplicated(coldata), ]
			object@list_of_data[[x]] <- data
			object@list_of_row_data[[x]] <- rowdata
			object@list_of_col_data[[x]] <- coldata
			return(object)
		})
		out <- add_function_info(object=out, function_name="mod_remove_duplicates", params=list(param="removes_duplicates"))
	 	return(out)
	})


