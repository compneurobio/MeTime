
#' Function to extract common samples across all datasets and store them only 
#' @description Modification(mod) function applied on object of class metime_analyser to extract common samples across datasets.
#' @examples # extracting common samples across all datasets
#' new_object_with_only_common_samples <- mod_extract_common_samples(object=metime_analyser_object)
#' @param object An object of class metime_anaylser
#' @return metime_analyser object with only common samples across all datasets present in the object parsed
#' @export
setGeneric("mod_extract_common_samples", function(object) standardGeneric("mod_extract_common_samples") )

setMethod("mod_extract_common_samples", "metime_analyser", function(object) {
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


