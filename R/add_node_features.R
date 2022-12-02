#' Function to add features to visnetwork plot from another plotter object
#' @description Function to add node features to see the nodes in the network that affected differently
#' @param network_plotter_object plotter object with network information
#' @param guide_plotter_object guide from which the colors are to be extracted
#' @param which_type type of the guide plotter object to be used. Current options are "regression" and "conservation"
#' @param metab_colname name of the column in guide plotter object that represents the metabolites
#' @return network plotter object with new node colors/features
#' @export
setGeneric("add_node_features", function(network_plotter_object, guide_plotter_object, which_type, metab_colname) standardGeneric("add_node_features"))
setMethod("add_node_features", "metime_plotter", function(network_plotter_object, guide_plotter_object, which_type, metab_colname) {
		network_plotter_object@plot_data[["node"]] <- network_plotter_object@plot_data[["node"]][order(network_plotter_object@plot_data[["node"]]$label), ]
		guide_plotter_object@plot_data[[1]] <- guide_plotter_object@plot_data[[1]][guide_plotter_object@plot_data[[1]][ ,metab_colname] %in% network_plotter_object@plot_data[[1]]$label, ]
		guide_plotter_object@plot_data[[1]] <- guide_plotter_object@plot_data[[1]][order(guide_plotter_object@plot_data[[1]][ ,metab_colname]), ]
		if(which_type %in% "regression") {
			column_for_colors <- guide_plotter_object@plot_data[[1]][ ,c("beta", "pval")]
			column_for_colors <- sign(column_for_colors$beta) * -log10(column_for_colors$pval)
		} else if(which_type %in% "conservation") {
			column_for_colors <- guide_plotter_object@plot_data[[1]][ ,"ci"]
		}
		color.gradient <- function(x, colors=c("blue","gray","red"), colsteps=50) {
  			return(colorRampPalette(colors) (colsteps)[findInterval(x, seq(-1, 1, length.out=colsteps))] )
		}
		gradient <- color.gradient(column_for_colors)
		network_plotter_object@plot_data[["node"]]$color <- gradient
		return(network_plotter_object)
	}) 
