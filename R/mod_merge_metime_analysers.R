
#' Function to merge one or more metime_analyser objects
#' @description function to merge multiple metime_analyser objects
#' @param list_of_objects list of metime analyser objects that are to be merged. The objects must be named
#' Ex: list_of_objects <- list(object1=object1, object2=object2, ...). If not named then they will be named by the function
#' @param annotations_index new list with annotations_index. Can also set to be NULL.
#' @returns A merged metime_analyser object
#' @export
mod_merge_metime_analysers <- function(list_of_objects, annotations_index) {
				if(is.null(names(list_of_objects))) {
					names(list_of_objects) <- paste("object_", 1:length(list_of_objects), sep="")
				}
				list_of_data <- lapply(seq_along(list_of_objects), function(index) {
						names(list_of_objects[[index]]@list_of_data) <- paste(names(list_of_objects)[index], "_", names(list_of_objects[[index]]@list_of_data))
						return(list_of_objects[[index]]@list_of_data)
					})
				list_of_col_data <- lapply(seq_along(list_of_objects), function(index) {
						names(list_of_objects[[index]]@list_of_col_data) <- paste(names(list_of_objects)[index], "_", names(list_of_objects[[index]]@list_of_col_data))
						return(list_of_objects[[index]]@list_of_col_data)
					})
				list_of_row_data <- lapply(seq_along(list_of_objects), function(index) {
						names(list_of_objects[[index]]@list_of_row_data) <- paste(names(list_of_objects)[index], "_", names(list_of_objects[[index]]@list_of_row_data))
						return(list_of_objects[[index]]@list_of_row_data)
					})
				list_of_results <- lapply(seq_along(list_of_objects), function(index) {
						names(list_of_objects[[index]]@results) <- paste(names(list_of_objects)[index], "_", names(list_of_objects[[index]]@results))
						return(list_of_objects[[index]]@results)
					})

				if(is.null(annotations_index)) {
					final_object <- new("metime_analyser", list_of_data=list_of_data, 
						list_of_row_data=list_of_row_data, list_of_col_data=list_of_col_data, annotations=NULL, 
						results=list_of_results)
				} else {
					final_object <- new("metime_analyser", list_of_data=list_of_data, 
						list_of_row_data=list_of_row_data, list_of_col_data=list_of_col_data, annotations=annotations_index,
						results=list_of_results)	
				}
				out <- final_object
				return(out)
	}

