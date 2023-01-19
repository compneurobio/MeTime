
#' Function to merge two sets of data for any analysis 
#' @description A function to merge two or more datasets and use it for analysis
#' @param object An S4 object of class metime_analyser
#' @param which_data Datasets to be merged. Only two or more are allowed
#' @param name character vector to define the name of the new dataset
#' @return A new S4 object of class metime_analyser with the new merged dataset appended to it
#' @export 
setGeneric("mod_merge_data", function(object, which_data, name) standardGeneric("mod_merge_data")) 
setMethod("mod_merge_data", "metime_analyser", function(object, which_data, name)	{
		stopifnot(length(which_data)>=2)
		new_object <- mod_extract_common_samples(object)
		test <- lapply(seq_along(which_data), function(i) {
				if(i==1) {
					data <- new_object@list_of_data[[which_data[i]]]
					row_data <- new_object@list_of_row_data[[which_data[i]]]
					col_data <- new_object@list_of_col_data[[which_data[i]]]
					data <- data[order(rownames(data)), ]
					row_data <- row_data[order(rownames(row_data)), ]
				} else {
					dummy_data <- new_object@list_of_data[[which_data[i]]]
					dummy_data <- dummy_data[order(rownames(dummy_data)), ]
					data <- cbind(data, dummy_data)
					dummy_row <- new_object@list_of_row_data[[which_data[i]]]
					dummy_row <- dummy_row[rownames(dummy_row) %in% rownames(row_data), ]
					dummy_row <- dummy_row[order(rownames(dummy_row)), ]
					dummy_row <- dummy_row %>% dplyr::select(-id, -subject, -time)
					row_data <- cbind(row_data, dummy_row)
					dummy_col <- object@list_of_col_data[[which_data[i]]]
					col_data <- plyr::rbind.fill(col_data, dummy_col)
				}
			})
		out <- get_append_analyser_object(object=object, data=data, col_data=col_data, row_data=row_data, name=name)
		out <- add_function_info(object=out, function_name="mod_merge_data", params=list(which_data=which_data, name=name))
		return(out)
})


