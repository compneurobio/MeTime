#' Function to get common samples at multiple timepoints chosen
#' @description Get a vector of common samples that can be used for stratification
#' @param object An S4 object of class metime_analyser
#' @param which_data Dataset/s of interest
#' @param timepoints A character vector to define timepoints
#' @importClassesFrom metime_analyser
#' @returns character vector of common subjects at common timepoints
#' @export
setGeneric("get_common_samples_at_timepoints", function(object, which_data, timepoints) standardGeneric("get_common_samples_at_timepoints"))
setMethod("get_common_samples_at_timepoints", "metime_analyser", function(object, which_data, timepoints) {
		if(length(which_data)==1) {
      		row_data <- get_rowdata(object, which_data=which_data)
		} else {
			object <- mod_extract_common_samples(object)
			object@list_of_data <- lapply(object@list_of_data, function(x) return(x[order(rownames(x)), ]))
			list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
			data <- do.call(cbind, unname(list_of_data))
			row_data <- object@list_of_row_data[[which_data[1]]]
			row_data <- row_data[rownames(row_data) %in% rownames(data), ]
		}
		row_data <- row_data[row_data$time %in% timepoints, ]
		row_data$subject %>% unique() %>% return()
	})