
#' Function to merge two sets of results
#' @description Modification(mod) function to merge two sets of results based on calc_type
#' @param object An S4 object of class metime_analyser
#' @param results_index character of numeric vector to define the index of results to merge.
#' if length of results_index is 1, then the plot_data list is merged into a single data.frame
#' @param sub_results Vector to define indices of plot_data to be merged. Set to 1 by default. 
#' @param groups character vector to define the two sets of different results that are being merged
#' @param name character to name the new merged results
#' @returns metime_analyser object with merged results 
#' @seealso [mod_merge_data]  
#' @export
mod_merge_results <- function(object, results_index, sub_results=1, groups, name) {
		if(length(results_index)==1) {
			results <- object@results[[results_index]]
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
				if(is.null(sub_results)) {
					warning("length of groups and results dataframes are not equal")
				} else {
					stopifnot(length(sub_results)==length(groups))
					new_data <- lapply(sub_results, function(b) {
							results$plot_data[[b]]$groups <- rep(groups[b], 
									each=length(results$plot_data[[b]][,1]))
							return(results$plot_data[[b]])	
						}) %>% do.call(what=rbind.data.frame)
					out <- get_make_results(object=object, data=list(new_data), metadata=NULL, 
						calc_info= paste(results$information$calc_info[sub_results], collapse=" || "),
						calc_type= paste(results$information$calc_type[sub_results[1]]),
						name=name)
					out <- add_function_info(object=out, function_name="mod_merge_results", 
						params=list(results_index=results_index, groups=groups, sub_results=sub_results))
					return(out)
				} 
			}
		} else {
			stopifnot(length(results_index)==length(groups))
			results_1 <- object@results[[results_index[1]]]
			results_2 <- object@results[[results_index[2]]]
			results_1_df <- results_1$plot_data[[1]] %>% as.data.frame()
			results_2_df <- results_2$plot_data[[1]] %>% as.data.frame()
			results_1_df$groups <- rep(groups[1], each=length(results_1_df[,1]))
			results_2_df$groups <- rep(groups[2], each=length(results_2_df[,1]))
			final_results <- rbind.data.frame(results_1_df, results_2_df)
			out <- get_make_results(object=object, data=list(new_data), metadata=NULL, 
						calc_info= paste(results_1$information$calc_info, results_2$information$calc_info, sep=" || "),
						calc_type= paste(results_1$information$calc_type),
						name=name)
			out <- add_function_info(object=out, function_name="mod_merge_results", 
						params=list(results_index=results_index, groups=groups, sub_results=sub_results))
			return(out)
		}
	}


