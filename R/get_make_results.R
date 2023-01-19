
#' Function to make a plottable object for viz functions
#' @description function to generate metime_plotter object from plot data and metadata
#' @param data dataframe of plotable data obtained from any calc object
#' @param metadata dataframe with the metadata for the plot table mentioned above. To obtain these see
#' get_metadata_for_rows() and get_metadata_for_columns()
#' @param calc_type A character to specify type of calculation - will be used for comp_ functions
#' For networks the accepted notations are "genenet_ggm", "multibipartite_ggm", and "temporal_network"
#' @param calc_info A string to define the information about calculation
#' @param plot_type type of the plot you want to build. eg: "box", "dot" etc. Its a character vector
#' @param style Style of plot, accepted inputs are "ggplot", "circos" and "visNetwork". Is a singular option.
#' @export
get_make_results <- function(object, data, metadata, calc_type, calc_info, plot_type, style, aesthetics) {
			plot_data <- list()
			empty_plots <- list()
			if(style %in% "visNetwork") {
				nodes <- unique(c(data$node1, data$node2))
    			node_list <- data.frame(id=1:length(nodes), label=nodes, group=as.character(1:length(nodes)))
    			for(i in 1:length(node_list$label)) {
          			g <- metadata[as.character(metadata$name) %in% as.character(node_list$label[i]), 2]
           			node_list$group[i] <- g
    			}
    			#Getting edge list
    			edge_list <- data.frame(from=1:length(data$node1), to=1:length(data$node2))
    			for(i in 1:length(data$node1)) {
        			edge_list$from[i] <- node_list[as.character(node_list$label) %in% as.character(data$node1[i]), "id"]
        			edge_list$to[i] <- node_list[as.character(node_list$label) %in% as.character(data$node2[i]), "id"]
   			 	}
   			 	if(calc_type %in% "genenet_ggm") {
   			 			dashes <- ifelse(data$pcor > 0, FALSE, TRUE)
   			 			edge_list$dashes <- dashes
    					edge_list$values <- data$pcor
   			 	} else if(calc_type %in% "multibipartite_ggm") {
   			 			dashes <- ifelse(data$coeffs.1 > 0 & data$coeffs.2 > 0, FALSE, TRUE)
   			 			edge_list$dashes <- dashes
    					edge_list$values <- (data$coeffs.1 + data$coeffs.2)/2 
    					edge_list$title <- paste("coeff1: ", data$coeffs.1, "<br /> coeff2: ", data$coeffs.2, sep=" ")
   			 	} else if(calc_type %in% "temporal_network") {
   			 			dashes <- ifelse(data$coeffs > 0, FALSE, TRUE)
   			 			edge_list$dashes <- dashes
   			 			edge_list$arrows <- rep("from", each=length(edge_list$dashes))
   			 	}
					object@results[[length(object@results)]][["plot_data"]][["node"]] <- node_list
					object@results[[length(object@results)]][["plot_data"]][["edge"]] <- edge_list
					object@results[[length(object@results)]][["plot_data"]][["metadata"]] <- metadata
			} else {
					if(is.null(metadata)) {
						plot_data[[1]] <- as.data.frame(data)
					} else {
						data <- data[order(rownames(data)), ]
						metadata <- metadata[rownames(metadata) %in% rownames(data), ]
						data <- data[rownames(data) %in% rownames(metadata),]
						object@results[[length(object@results)]][["plot_data"]] <- as.data.frame(cbind(data, metadata))
					}
			}
			if(style %in% "ggplot") {
				if(plot_type %in% "dot") {
					object@results[[length(object@results)]][["plots"]] <- ggplot(object@results[[length(object@results)]][["plot_data"]], aes_string(x=aesthetics$x, y=aesthetics$y)) +
												geom_point()
				} else if(plot_type %in% "heatmap") {
					object@results[[length(object@results)]][["plot"]] <- ggplot(object@results[[length(object@results)]][["plot_data"]], aes_string(x=aesthetics$x, y=aesthetics$y)) +
												geom_tile(aes_string(fill=aesthetics$fill))
				} else {
					object@results[[length(object@results)]][["plot"]] <- ggplot(object@results[[length(object@results)]][["plot_data"]])
				}
			} else if(style %in% "circos") {
						empty_plots[[i]] <- NULL
			} else if(style %in% "visNetwork") {
						empty_plots[[i]] <- NULL
			}
			return(object)
}

