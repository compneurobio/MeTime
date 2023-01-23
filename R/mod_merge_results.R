
#' Function to combine plotter objects
#' @description Function to combine plotter objects based on similar calc type
#' @param object An S4 object of class metime_analyser
#' @param results_indices one or more indices of results to merge data
#' @param groups character vector to define the groups that are involved
#' @return object with add merged letters  
#' @export
mod_merge_results <- function(object, results_indices, groups) {
			results <- object@results[results_indices]
			if(length(results) > 1) {
				types <- vapply(seq_along(results), function(x) {
						type <- results[[x]]$information$calc_type
						return(type)
					}, character(1))
				if(length(unique(types))==1) {
					data <- lapply(seq_along(results), function(x) {
							new_data <- lapply(seq_along(results[[x]]$plot_data), function(y) {
									dummy_data <- results[[x]]$plot_data[[y]]
									dummy_data$group <- rep(groups[y], each=length(dummy_data[,1]))
									return(dummy_data)
								})
						})
				} else {
					stop("Results are not of the same type")
				}
			} else {
				if(length(groups)!=length(results$plot_data)) stop("number of groups are not equal")
			}
	}

