
#' Function to remove NA's from data matrices
#' @description Modification(mod) function applied on S4 object to remove NA's and change col_data and row_data accordingly
#' @param object S4 object of class metime_analyser
#' @param which_data character vector to define the datasets to be used. Can also use a numeric vector with indices.
#' @return S4 object with NA's removed and col_data, row_data manipulated accordingly
#' @export
setGeneric("mod_remove_nas", function(object, which_data) standardGeneric("mod_remove_nas"))
setMethod("mod_remove_nas", "metime_analyser", function(object, which_data) {
				if(!all(which_data %in% names(object@list_of_data))) {
					warning("The datasets mentioned are not found in the metime_analyser object. Exiting without making any changes.")
					return(object)
				}
				object@list_of_data[names(object@list_of_data) %in% which_data] <- lapply(object@list_of_data[names(object@list_of_data) %in% which_data], 
																	function(x) return(na.omit(x)))
				rownames <- lapply(seq_along(object@list_of_data[names(object@list_of_data) %in% which_data]), function(a) {
						object@list_of_data[names(object@list_of_data) %in% which_data][[a]] %>% rownames() %>% return()
					})
				colnames <- lapply(seq_along(object@list_of_data[names(object@list_of_data) %in% which_data]), function(a) {
						object@list_of_data[names(object@list_of_data) %in% which_data][[a]] %>% colnames() %>% return()
					})
				object@list_of_row_data[names(object@list_of_row_data) %in% which_data] <- lapply(seq_along(object@list_of_row_data[names(object@list_of_row_data) %in% which_data]), function(x) {
						rows <- object@list_of_row_data[names(object@list_of_row_data) %in% which_data][[x]]
						rows <- rows[rownames(rows) %in% rownames[[x]], ]
						return(rows)
					})
				object@list_of_col_data[names(object@list_of_col_data) %in% which_data] <- lapply(seq_along(object@list_of_col_data[names(object@list_of_col_data) %in% which_data]), function(x) {
						cols <- object@list_of_col_data[names(object@list_of_col_data) %in% which_data][[x]]
						cols <- cols[rownames(cols) %in% colnames[[x]], ]
						return(cols)
					})
				names(object@list_of_col_data) <- names(object@list_of_data)
				out <- object
				out <- add_function_info(object=out, function_name="mod_remove_nas", params=list(which_data=which_data))
				return(out)
	})


