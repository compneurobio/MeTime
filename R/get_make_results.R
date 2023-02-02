
#' Function to make results list for metime_analyser object
#' @description function to generate results for metime_analyser object
#' @param object An S4 object of class metime_analyser
#' @param data list of dataframes of plotable data obtained from any calc function
#' @param metadata dataframe or a list of dataframes with the metadata for the plot table mentioned above. To obtain these see
#' get_metadata_for_rows() and get_metadata_for_columns()
#' @param calc_type A character vector to specify type of calculation - will be used for comp_ functions
#' For networks the accepted notations are "genenet_ggm", "multibipartite_ggm", and "temporal_network"
#' sjould be the same length as the list of data provided
#' @param calc_info A string to define the information about calculation, should be the same length as the list
#' data provided
#' @param name Name of the result 
#' @return object with results of the calculation updated
#' @export
setGeneric("get_make_results", function(object, data, metadata, calc_type, calc_info, name) standardGeneric("get_make_results"))
setMethod("get_make_results", "metime_analyser", function(object, data, metadata, calc_type, calc_info, name) {
			stopifnot(length(calc_type)==length(calc_info))
			stopifnot(length(calc_type)==length(data))
			stopifnot(length(calc_info)==length(data))
			if(length(grep("ggm|network", calc_type))==1) {
				data <- data[[1]]
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
   			 	if(length(grep("calc_|mod_merge_results", names(object@results[[length(object@results)]]$functions_applied))) ==1) {
					object@results[[length(object@results)+1]] <- list(functions_applied=list(), 
						plot_data=list(node=node_list, edge=edge_list, metadata=metadata),
						information=list(calc_type=calc_type, calc_info=calc_info))
				} else {
					object@results[[length(object@results)]]$plot_data$node <- node_list
					object@results[[length(object@results)]]$plot_data$edge <- edge_list
					object@results[[length(object@results)]]$plot_data$metadata <- metadata
					object@results[[length(object@results)]]$information$calc_type <- calc_type
					object@results[[length(object@results)]]$information$calc_info <- calc_info
				}
			} else {
				if(is.null(metadata)) {
					plot_data <- data
				} else {
					plot_data <- lapply(seq_along(data), function(x) {
							if(length(metadata)==0) {
								return(data[[x]])
							}
							if(class(metadata) %in% "list") {	
								data[[x]] <- data[[x]][order(rownames(data[[x]])), ]
								dummy_metadata <- metadata[[x]][rownames(metadata[[x]]) %in% rownames(data[[x]]), ]
								data[[x]] <- data[[x]][rownames(data[[x]]) %in% rownames(dummy_metadata), ]
							} else {
								data[[x]] <- data[[x]][order(rownames(data[[x]])), ]
								dummy_metadata <- metadata[rownames(metadata) %in% rownames(data[[x]]), ]
								data[[x]] <- data[[x]][rownames(data[[x]]) %in% rownames(dummy_metadata), ]
							}
							return(cbind.data.frame(data[[x]], dummy_metadata))
					})
				}
				if(length(grep("calc_|mod_merge_results", names(object@results[[length(object@results)]]$functions_applied)))==1) {
					object@results[[length(object@results)+1]] <- list(functions_applied=list(), plot_data=plot_data,
											information=list(calc_type=calc_type, calc_info=calc_info))
				} else {
					object@results[[length(object@results)]]$plot_data <- plot_data
					object@results[[length(object@results)]]$information$calc_type <- calc_type
					object@results[[length(object@results)]]$information$calc_info <- calc_info
				}
			}
			names(object@results)[length(object@results)] <- name
			return(object)
})


