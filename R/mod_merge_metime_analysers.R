
#' Function to merge one or more metime_analyser objects
#' @description function to merge multiple metime_analyser objects. Will not be displayed in add_function_info()
#' @param ... metime analyser objects that are to be merged. The objects can be named by parsing it as an argument
#' If not named then they will be named by the function
#' @param annotations_index new list with annotations_index. If set to NULL the first object's annotations are taken with modification.
#' @seealso [add_function_info]
#' @returns A merged metime_analyser object
#' @export
mod_merge_metime_analysers <- function(..., annotations_index) {
				list_of_objects <- list(...)
				if(is.null(names(list_of_objects))) {
					names(list_of_objects) <- paste("object_", 1:length(list_of_objects), sep="")
				}
				list_of_data <- lapply(seq_along(list_of_objects), function(index) {
						names(list_of_objects[[index]]@list_of_data) <- paste(names(list_of_objects)[index], "_", names(list_of_objects[[index]]@list_of_data), sep="")
						return(list_of_objects[[index]]@list_of_data)
					})
				list_of_col_data <- lapply(seq_along(list_of_objects), function(index) {
						names(list_of_objects[[index]]@list_of_col_data) <- paste(names(list_of_objects)[index], "_", names(list_of_objects[[index]]@list_of_col_data), sep="")
						return(list_of_objects[[index]]@list_of_col_data)
					})
				list_of_row_data <- lapply(seq_along(list_of_objects), function(index) {
						names(list_of_objects[[index]]@list_of_row_data) <- paste(names(list_of_objects)[index], "_", names(list_of_objects[[index]]@list_of_row_data), sep="")
						return(list_of_objects[[index]]@list_of_row_data)
					})
				list_of_results <- lapply(seq_along(list_of_objects), function(index) {
						names(list_of_objects[[index]]@results) <- paste(names(list_of_objects)[index], "_", names(list_of_objects[[index]]@results), sep="")
						return(list_of_objects[[index]]@results)
					})

				if(is.null(annotations_index)) {
					annotations_index <- list_of_objects[[1]]@annotations
					annotations_index <- lapply(seq_along(annotations_index), function(x) {
							annotations_index[[x]] <- paste(names(list_of_objects[1]), "_", annotations_index[[x]], sep="")
						})
					final_object <- new("metime_analyser", list_of_data=list_of_data, 
						list_of_row_data=list_of_row_data, list_of_col_data=list_of_col_data, annotations=annotations_index, 
						results=list_of_results)
				} else {
					final_object <- new("metime_analyser", list_of_data=list_of_data, 
						list_of_row_data=list_of_row_data, list_of_col_data=list_of_col_data, annotations=annotations_index,
						results=list_of_results)	
				}
				out <- final_object
				return(out)
	}

