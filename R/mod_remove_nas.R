
#' Function to remove NA's from data matrices
#' @description Modification(mod) function applied on S4 object to remove NA's and change col_data and row_data accordingly
#' @param object S4 object of class metime_analyser
#' @param which_data character vector to define the datasets to be used. Can also use a numeric vector with indices.
#' @return S4 object with NA's removed and col_data, row_data manipulated accordingly
#' @export
setGeneric("mod_remove_nas", function(object, which_data) standardGeneric("mod_remove_nas"))
setMethod("mod_remove_nas", "metime_analyser", function(object, which_data) {
				list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
				list_of_data <- lapply(list_of_data, function(x) return(na.omit(x)))
				object@list_of_row_data <- lapply(seq_along(object@list_of_row_data), function(x) {
						if(names(object@list_of_row_data)[x] %in% which_data) {
							rows <- object@list_of_row_data[[x]]
							rows <- rows[rownames(rows) %in% rownames(list_of_data[[x]]), ]
							return(rows)
						} else {
							return(object@list_of_row_data[[x]])
						}
					})
				names(object@list_of_row_data) <- names(object@list_of_data)
				object@list_of_col_data <- lapply(seq_along(object@list_of_col_data), function(x) {
						if(names(object@list_of_col_data)[x] %in% which_data) {
							cols <- object@list_of_col_data[[x]]
							cols <- cols[rownames(cols) %in% colnames(list_of_data[[x]]), ]
							return(cols)
						} else {
							return(object@list_of_col_data[[x]])
						}
					})
				names(object@list_of_col_data) <- names(object@list_of_data)
				out <- object
				out <- add_function_info(object=out, function_name="mod_remove_nas", params=list(which_data=which_data))
				return(out)
	})


