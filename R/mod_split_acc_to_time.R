

#' Function to split data acoording to time
#' @description Function to split the list of dataframes into a nested list with each dataframe 
#' being split into into dataframes of different timepoints
#' @examples #splitting data according to time
#' new_data <- mod_split_acc_to_time(object=metime_analyser_object)
#' @param object An object of class metime_analyser
#' @return list_of_data with each dataframe being broken into a list of dataframes with respect to the timepoint they belong to
#' @export
setGeneric("mod_split_acc_to_time", function(object) standardGeneric("mod_split_acc_to_time") )

setMethod("mod_split_acc_to_time", "metime_analyser", function(object) {
	list_of_data <- object@list_of_data
	list_of_data_temporals <- lapply(list_of_data, function(data) {
		names <- strsplit(rownames(data), split="_")
		times <- unlist(lapply(names, function(x) return(x[2])))
		indList <- split(seq_along(times), times)

		#timepoint separation of the data

		list_of_temporals <- lapply(indList, function(x) {
				return(data[x,])
		})
		names(list_of_temporals) <- names(indList)
		return(list_of_temporals)
	})
	names(list_of_data_temporals) <- names(list_of_data)
	object@list_of_data <- list_of_data_temporals
	out <- object
	return(out)
})


