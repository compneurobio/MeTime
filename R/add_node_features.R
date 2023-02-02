#' Function to add features to visnetwork plot from another plotter object
#' @description Function to add node features to see the nodes in the network that affected differently
#' @param object An S4 object of class metime_analyser
#' @param results_indices indices as a list to define which results to use. Eg: list(network=1/"name of the results", guide=2)
#' @param which_calculation index for plot_data to be used. Set to 1 by default.
#' @param metab_colname name of the column in guide plotter object that represents the metabolites. set to "name" by default
#' @return network plotter object with new node colors/features
#' @export
setGeneric("add_node_features", function(object, results_indices, which_calculation=1, metab_colname="name") standardGeneric("add_node_features"))
setMethod("add_node_features", "metime_plotter", function(object, results_indices, which_calculation=1, metab_colname="name") {
		stopifnot(all(names(results_indices) %in% c("network", "guide")))
		stopifnot(is.null(results_indices$network))
		stopifnot(is.null(results_indices$guide))
		color.gradient <- function(x, colors=c("blue","gray","red"), colsteps=50) {
  				return(colorRampPalette(colors) (colsteps)[findInterval(x, seq(-1, 1, length.out=colsteps))] )
		}
		network <- object@results[[results_indices$network]]
		guide <- object@results[[results_indices$guide]]
		data_of_interest <- guide$plot_data[[which_results]]
		data_of_interest <- data_of_interest[order(data_of_interest[ ,metab_colname]), ]
		network$plot_data$node <- network$plot_data$node[order(network$plot_data$node$label), ]
		data_of_interest <- data_of_interest[data_of_interest[ ,metab_colname] %in% network$node$label, ]
		if(guide$information$calc_type[which_results] %in% "regression") {
			column_for_colors <- data_of_interest[ ,c("beta", "pval")]
			column_for_colors <- sign(column_for_colors$beta) * -log10(column_for_colors$pval)
		} else if(guide$information$calc_type[which_results] %in% "CI_metabolite") {
			column_for_colors <- data_of_interest[ ,"ci"]
		}
		gradient <- color.gradient(column_for_colors)
		network$plot_data$node$color <- gradient
		object@results[[results_indices$network]] <- network
		return(object)
	}) 


