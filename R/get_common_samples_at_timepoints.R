#' Get common samples at multiple timepoints chosen
#' @description Get a vector of common samples that can be used for stratification.
#' @param objecta S4 object of the class "metime_analyzer".
#' @param which_data a character to define which dataset is to be used.
#' @param timepoints a character vector to define time points
#' @returns A character vector of common subjects at common time points.
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
		row_data_list <- lapply(timepoints, function(tp) {
				row_data[row_data$time %in% tp, ] %>% return()
			})
		list_of_samples <- lapply(seq_along(row_data_list), function(ind) row_data_list[[ind]]$subject %>% return())
		common_samples <- Reduce(intersect, list_of_samples) 
		row_data <- row_data %>% filter(subject %in% common_samples)
		row_data$subject %>% unique() %>% return()
	})