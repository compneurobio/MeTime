#' Function to remove datasets or results from the metime analyser object
#' @description An S4 method to remove unwanted data
#' @param object An S4 object of class metime_analyser
#' @param which_data index or name of the dataset to be removed or result
#' @param type "result" or "dataset" based on the type to be removed
#' @returns object with results/dataset removed
#' @export
setGeneric("rm_data", function(object, which_data, type="dataset") standardGeneric("rm_data"))
setMethod("rm_data", "metime_analyser", function(object, which_data, type="dataset") {
		if(length(which_data)==0) {
			warning("which_data is not specified. Exiting without making any changes")
			return(object)
		}
		if(type %in% "dataset") {
			object@list_of_data[[which_data]] <- NULL
			object@list_of_col_data[[which_data]] <- NULL
			object@list_of_row_data[[which_data]] <- NULL
		} else if(type %in% "results") {
			object@results[[which_data]] <- NULL
		} else {
			warning("Please check the input for type. Exiting without making any changes")
		}
		return(object)
	})