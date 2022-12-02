
#' Function to get only common samples from the dataframes in list_of_data 
#' @description A method applied on object of class metime_analyse to extract common samples across datasets. Also has an option to split the data according 
#' timepoints(mod_split_acc_time()).
#' @examples # extracting common samples across all datasets
#' new_list_of_data <- mod_common_sample_extractor(object=metime_analyser_object)
#' @param object An object of class metime_anaylser
#' @param time_splitter A boolean input: True leads to splitting of the data wrt time, 
#'	 					False returns all the dataframes as they are with common rows
#' @return list_of_data with common samples across all time points
#' @export
setGeneric("mod_extract_common_samples", function(object, time_splitter=FALSE) standardGeneric("mod_extract_common_samples") )

setMethod("mod_extract_common_samples", "metime_analyser",function(object, time_splitter=FALSE) {
		list_of_names <- lapply(object@list_of_data, function(x) {
					return(rownames(x))
			})
		common_samples <- Reduce(intersect, list_of_names)
		object@list_of_data <- lapply(object@list_of_data, function(x) {
					x <- x[rownames(x) %in% common_samples, ]
					x <- x[order(rownames(x)), ]
					return(x)
			})
		if(time_splitter) {
				object <- mod_split_acc_to_time(object)
		}
		out <- object
		return(out)
})


