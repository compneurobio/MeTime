
#' Function to extract col data of a dataset
#' @description Function to get coldata 
#' @param object An object of class S4 
#' @param which_data Dataset of interest
#' @return col data of the dataset of interest
#' @export
setGeneric("get_coldata", function(object, which_data) standardGeneric("get_coldata"))
setMethod("get_coldata", "metime_analyser", function(object, which_data) {
				data <- object@list_of_col_data[[which_data]]
				return(as.data.frame(data))
	})