
#' Function to merge one or more metime_analyser objects
#' @description function to merge multiple metime_analyser objects
#' @param list_of_objects list of metime analyser objects that are to be merged
#' @param annotations_index new list with annotations_index. Can also set to be NULL.
#' @returns A merged metime_analyser object
#' @export
setGeneric("mod_merge_metime_analysers", function(list_of_objects, annotations_index) standardGeneric("mod_merge_metime_analysers"))
setMethod("mod_merge_metime_analysers", "metime_analyser", function(list_of_objects, annotations_index) {
				list_of_data <- lapply(list_of_objects, function(x) return(x@list_of_data))
				list_of_col_data <- lapply(list_of_objects, function(x) return(x@list_of_col_data))
				list_of_row_data <- lapply(list_of_objects, function(x) return(x@list_of_row_data))
				if(is.null(annotations_index)) {
						final_object <- new("metime_analyser", list_of_data=list_of_data, 
							list_of_row_data=list_of_row_data, list_of_col_data=list_of_col_data, annotations=NULL)
				} else {
						final_object <- new("metime_analyser", list_of_data=list_of_data, 
							list_of_row_data=list_of_row_data, list_of_col_data=list_of_col_data, annotations=annotations_index)	
				}
				out <- final_object
				return(out)
	}) 

