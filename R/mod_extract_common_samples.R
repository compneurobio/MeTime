
#' Function to get only common samples from the dataframes in list_of_data 
#' @description A method applied on object of class metime_analyse to extract common samples across datasets.
#' @examples # extracting common samples across all datasets
#' new_list_of_data <- mod_extract_common_samples(object=metime_analyser_object)
#' @param object An object of class metime_anaylser
#' @return list_of_data with common samples across all time points
#' @export
setGeneric("mod_extract_common_samples", function(object) standardGeneric("mod_extract_common_samples") )

setMethod("mod_extract_common_samples", "metime_analyser",function(object) {
		list_of_names <- lapply(object@list_of_data, function(x) {
					return(rownames(x))
			})
		common_samples <- Reduce(intersect, list_of_names)
		object@list_of_data <- lapply(object@list_of_data, function(x) {
					x <- x[rownames(x) %in% common_samples, ]
					x <- x[order(rownames(x)), ]
					return(x)
			})
		object@list_of_row_data <- lapply(object@list_of_row_data, function(x) {
					x <- x[rownames(x) %in% common_samples, ]
					x <- x[order(rownames(x)), ]
					return(x)
			})
		out <- object
		out <- add_function_info(object=out, function_name="mod_extract_common_samples", 
			params=list(param1="no_additional_information"))
		return(out)
})


