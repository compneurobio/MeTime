
#' Function to merge one or more metime_analyser objects
#' @description function to merge multiple metime_analyser objects
#' @param list_of_objects list of metime analyser objects that are to be merged
#' @param annotations_index new list with annotations_index. Can also set to be NULL.
#' @importClassesFrom metime_analyser
#' @returns A merged metime_analyser object
#' @export
mod_merge_metime_analysers <- function(list_of_objects, annotations_index) {
				list_of_data <- lapply(list_of_objects, function(x) return(x@list_of_data))
				list_of_col_data <- lapply(list_of_objects, function(x) return(x@list_of_col_data))
				list_of_row_data <- lapply(list_of_objects, function(x) return(x@list_of_row_data))
				list_of_results <- lapply(list_of_objects, function(x) return(x@results))
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

