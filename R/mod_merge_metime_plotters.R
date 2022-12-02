
#' Function to combine plotter objects
#' @description Function to combine plotter objects based on similar calc type
#' @param list_of_objects list of plotter objects to be merged
#' @param groups character vector to define the groups that are involved
#' @return plotter object merged with all the plotters 
#' @export
mod_merge_metime_plotters <- function(list_of_objects, groups) {
			plotter_objects <- list_of_objects
			calc_types <- unlist(lapply(plotter_objects, function(x) return(x@calc_type)))
			if(length(unique(calc_types)) > 1) stop("Calc_types of all the plotter objects is not same")
			if(!(length(groups) == length(plotter_objects))) stop("length of plotter objects and groups don't match check the input")
			plotter_objects <- lapply(1:length(groups), function(x) {
						plotter_objects[[x]]@plot_data[[1]]$group <- rep(groups[x], each=length(plotter_objects[[x]]@plot_data[[1]][,1]))
						return(plotter_objects[[x]])
				})
			combined_plotter_data <- lapply(plotter_objects, function(x) return(x@plot_data[[1]])) %>% do.call(what=rbind.data.frame)
			empty_plot <- ggplot(data=combined_plotter_data)
			calc_type <- unique(calc_types)
			calc_info <- "merged_calculations"
			style <- plotter_objects[[1]]@style
			plot_type <- plotter_objects[[1]]@plot_type
			out <- new("metime_plotter", plot_data=list(combined_plotter_data), plot=list(empty_plot), calc_type=calc_type, 
						calc_info=calc_info, plot_type=plot_type, style=style)
			return(out)
	}

