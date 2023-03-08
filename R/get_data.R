#' Function to get particular dataset's data 
#' @description function to get data from S4 object
#' @param object An S4 object of class metime_analyser
#' @param which_data character/index to define the data of interest
#' @importClassesFrom metime_analyser
#' @returns data.frame with data of interest
#' @export
setGeneric("get_data", function(object, which_data) standardGeneric("get_data"))
setMethod("get_data", "metime_analyser", function(object, which_data) {
			stopifnot(!is.null(which_data))
			return(object@list_of_data[[which_data]] %>% as.data.frame())
	})