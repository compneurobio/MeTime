
#' Function to extract row data of a dataset
#' @description Function to get rowdata 
#' @param object An object of class S4 
#' @param which_data Dataset of interest
#' @return row data of the dataset of interest
#' @export
setGeneric("get_rowdata", function(object, which_data) standardGeneric("get_rowdata"))
setMethod("get_rowdata", "metime_analyser", function(object, which_data) {
				data <- object@list_of_row_data[[which_data]]
				return(as.data.frame(data))
	})


