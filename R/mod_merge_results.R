
#' Function to merge two sets of results
#' @description Modification(mod) function to merge two sets of results based on calc_type
#' @param object An S4 object of class metime_analyser
#' @param results_index character of numeric vector to define the index of results to merge.
#' if length of results_index is 1, then the plot_data list is merged into a single data.frame
#' @param sub_results List of character/numeric vectors to define indices of plot_data to be merged. 
#' This is not needed when trying to merge plot_data of same results but need it for merging different
#' results. The length of sub_results should be equal to length of results_index and also should be in
#' the same order. 
#' @param groups character vector to name the results you want to merge and should be
#' of length plot_data in the case of merging results of same calculation. 
#' Else groups should be of length equal to unlist(sub_results) and the groups should follow the order with
#' precedence to results_index.
#' Eg: results_index=1:2, sub_results=list(c(3,4), c(5,6)), then groups <- c("name_1_3", "name_1_4", "name_2_5", "name_2_6")
#' @param name character to name the new merged results
#' @returns metime_analyser object with merged results 
#' @seealso [mod_merge_data], [mod_merge_row_data_and_data]  
#' @export
mod_merge_results <- function(object, results_index, sub_results=1, groups, name) {
		if(length(results_index)==0) {
			warning("results_index needs to be of length >=1. Exiting without making any changes.")
			return(object)
		}
		if(length(results_index)==1) {
			results <- object@results[[results_index]]
			if(!is.null(results$plot_data$network)) {
				warning("Cannot merge results of type networks like this. Exiting without making any changes.")
				return(object)
			}
			if(length(results$plot_data) == length(groups)) {
				new_data <- lapply(seq_along(results$plot_data), function(b) {
						results$plot_data[[b]]$groups <- rep(groups[b], 
									each=length(results$plot_data[[b]][,1]))
							return(results$plot_data[[b]])	
					}) %>% do.call(what=rbind.data.frame)
				out <- get_make_results(object=object, data=list(new_data), metadata=NULL, 
						calc_info= paste(results$information$calc_info, collapse=" || "),
						calc_type= paste(results$information$calc_type[1]),
						name=name)
				out <- add_function_info(object=out, function_name="mod_merge_results", 
					params=list(results_index=results_index, groups=groups))
				return(out)
			} else {
				if(length(unlist(sub_results))!=length(groups)) {
					warning("length of groups and results dataframes are not equal. Exiting without making any changes.")
					return(object)
				} else {
					sub_indices <- unname(unlist(sub_results))
					new_data <- lapply(seq_along(sub_indices), function(ind) {
							b <- sub_indices[ind]
							results$plot_data[[b]]$groups <- rep(groups[ind], 
									each=length(results$plot_data[[b]][,1]))
							return(results$plot_data[[b]])	
						}) %>% do.call(what=rbind.data.frame)
					out <- get_make_results(object=object, data=list(new_data), metadata=NULL, 
						calc_info= paste(results$information$calc_info[sub_indices], collapse=" || "),
						calc_type= paste(results$information$calc_type[sub_indices[1]]),
						name=name)
					out <- add_function_info(object=out, function_name="mod_merge_results", 
						params=list(results_index=results_index, groups=groups, sub_results=sub_results))
					return(out)
				} 
			}
		} else {
			if(length(results_index)!=length(sub_results)) {
				warning("Length of sub_results list and results_index is not the same. Are you trying to merge two different results or same results? Exiting without making any changes.")
				return(object)
			}
			if(length(unlist(sub_results))!=length(groups)) {
				warning("Length of groups don't match the length of results you want to merge. Exiting without making any changes.")
				return(object)
			}
			results_list <- lapply(results_index, function(ind) {
					object@results[[ind]]$plot_data
				})
			calc_types <- lapply(results_index, function(x) {
					object@results[[x]]$information$calc_type
				}) %>% unlist(recursive=FALSE) %>% unique()
			if(length(calc_types)!=1) {
				warning("Please check the type of calculations you are merging. Exiting without making any changes.")
				return(object)
			}
			if(!(grep("ggm|network", calc_types) %>% length()==1)) {
				results_combined <- lapply(seq_along(results_list), function(ind) {
						results_list[[ind]][sub_results[[ind]]]
					}) %>% unlist(recursive=FALSE)
			} else {
				results_combined <- lapply(seq_along(results_list), function(ind) {
						results_list[[ind]]$network$edge
					})
				calc_types <- paste("merged_", calc_types, sep="")
			}
			results_combined <- lapply(seq_along(results_combined), function(ind) {
						results_combined[[ind]]$groups <- rep(groups[ind], each=length(results_combined[[ind]][,1]))
						return(results_combined[[ind]])
				}) %>% do.call(what=plyr::rbind.fill)
			
			out <- get_make_results(object=object, data=list(results_combined), metadata=NULL, 
						calc_info= "Merged results of different kinds. See functions_applied for better information.",
						calc_type= calc_types,
						name=name)
			out <- add_function_info(object=out, function_name="mod_merge_results", 
						params=list(results_index=results_index, groups=groups, sub_results=sub_results))
			return(out)
		}
	}

