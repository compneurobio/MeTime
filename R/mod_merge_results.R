
#' Function to combine plotter objects
#' @description Function to combine plotter objects based on similar calc type
#' @param object An S4 object of class metime_analyser
#' @param results_index index of results to merge results of similar type
#' @param sub_results Vector to define indices of plot_data to be merged. Set to NULL by default.
#' @param groups character vector to define the groups that are involved
#' @param name character to name the new merged results
#' @returns metime_analyser object with merged results   
#' @export
mod_merge_results <- function(object, results_index, sub_results=NULL, groups, name) {
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
					stop("length of groups and results dataframes are not equal. Please fill in the right values for sub_results")
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
	}


