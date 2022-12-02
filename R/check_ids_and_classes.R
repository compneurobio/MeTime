

#' Function to check the ids in the data and data format
#' @description sanity check to check for ids and order of the data
#' @param object S4 object of class of metime_analyser
#' @return NULL if it passes all the sanity checks
#' @export
setGeneric("check_ids_and_classes",  function(object) standardGeneric("check_ids_and_classes")) 
setMethod("check_ids_and_classes", "metime_analyser", function(object) {
			list_of_interest <- list(data=object@list_of_data, col_data=object@list_of_col_data, row_data=object@list_of_row_data)
			names <- names(object@list_of_data)
			type <- unlist(lapply(strsplit(names(object), split="list_of_"), function(x) return(x[2])))
			for(i in 1:length(list_of_interest)) {
				for(j in 1:length(list_of_interest[[i]])) {
					if(class(!list_of_interest[[i]][[j]]) %in% "data.frame") {
						stop(paste("In dataset", names[j], "datatype", type[i], 
							"class of data should be data.frame only"))
					} 
					if(!all(rownames(list_of_interest[[i]][[j]] == list_of_interest[[i]][[j]]$id))) {
						stop(paste("In dataset", names[j], "datatype", type[i], 
							"rownames and ids are not mapped correctly"))
					}
				}
			}
			return(NULL)
	})


