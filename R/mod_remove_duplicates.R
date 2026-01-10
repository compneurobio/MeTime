
#' Function to remove duplicates
#' @description Modification(mod) function to remove duplicates from the datasets
#' @param object An S4 object of class metime_analyser
#' @return object after removing duplicated data
#' @export
setGeneric("mod_remove_duplicates", function(object) standardGeneric("mod_remove_duplicates"))
setMethod("mod_remove_duplicates", "metime_analyser", function(object) {
		out <- object
		for(x in names(out@list_of_data)) {
			data <- out@list_of_data[[x]]
			rowdata <- out@list_of_row_data[[x]]
			coldata <- out@list_of_col_data[[x]]
			data <- data[!duplicated(data), ]
			rowdata <- rowdata[!duplicated(rowdata), ]
			coldata <- coldata[!duplicated(coldata), ]
			out@list_of_data[[x]] <- data
			out@list_of_row_data[[x]] <- rowdata
			out@list_of_col_data[[x]] <- coldata
		}
		out <- add_function_info(object=out, function_name="mod_remove_duplicates", params=list(param="removes_duplicates"))
	 	return(out)
	})

