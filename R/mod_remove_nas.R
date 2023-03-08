
#' Function to remove NA's from data matrices
#' @description A method applied on S4 object to remove NA's and change data accordingly
#' @param object S4 object of class metime_analyser
#' @param which_data dataset/s for which the method is to be applied
#' @importClassesFrom metime_analyser
#' @return S4 object with NA's removed and data manipulated accordingly
#' @export
setGeneric("mod_remove_nas", function(object, which_data) standardGeneric("mod_remove_nas"))
setMethod("mod_remove_nas", "metime_analyser", function(object, which_data) {
				list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
				list_of_data <- lapply(list_of_data, function(x) return(na.omit(x)))
				test <- lapply(seq_along(list_of_data), function(i) {
						rows <- object@list_of_row_data[[which_data[i]]]
						cols <- object@list_of_col_data[[which_data[i]]]
						rows <- rows[rownames(rows) %in% rownames(list_of_data[[i]]), ]
						cols <- cols[rownames(cols) %in% colnames(list_of_data[[i]]), ]
						object@list_of_row_data[[i]] <- rows
						object@list_of_col_data[[i]] <- cols					
				})
				out <- object
				out <- add_function_info(object=out, function_name="mod_remove_nas", params=list(which_data=which_data))
				return(out)
	})


