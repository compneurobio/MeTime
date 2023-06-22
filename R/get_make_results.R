#' Get results list for S4 object of class "metime_analyser" object
#' @description Compile a result element for metime_analyser object.
#' @param object a S4 object of the class "metime_analyzer".
#' @param data a list of dataframes of plotable data obtained from any calc function.
#' @param metadata a dataframe or a list of dataframes with the metadata for the plot table mentioned above. To obtain these see
#' get_metadata_for_rows() and get_metadata_for_columns().
#' @param calc_type a character vector to specify type of calculation. Should be the same length as the list data provided.
#' @param calc_info a character to define the information about calculation, should be the same length as the list data provided.
#' @param name a character to be used as the index of the result.
#' @return A S4 object of class "metime_analyser" with the result appended to the list of results.
#' @export
setGeneric("get_make_results", function(object, data, metadata, calc_type, calc_info, name) standardGeneric("get_make_results"))
setMethod("get_make_results", "metime_analyser", function(object, data, metadata, calc_type, calc_info, name) {
			if(all(length(calc_type)!=length(calc_info), length(calc_type)!=length(data), length(calc_info)!=length(data))) {
				warning("Length of calc_type, calc_info and plot_data is not equal. Exiting without making any changes")
				return(object)
			}
			if(length(grep("ggm|network", calc_type))==1 & length(grep("merged", calc_type))!=1) {
				#Getting data
				data <- data[[1]]
				#Getting node list along with metadata
				nodes_names <- unique(c(data$node1, data$node2))
    			nodes <- data.frame(id=1:length(nodes_names), label=nodes_names)
				nodes <- dplyr::left_join(nodes, metadata, by=c("label"="id"))
				nodes$title <- nodes$label
    			#Getting edge list
    			edges <- data %>%
  						dplyr::mutate(from = match(node1, nodes$label),
         						to = match(node2, nodes$label))
  				#Adjusting edge list based on the type of the network
   			 	if(calc_type %in% "genenet_ggm") {
   			 		dashes <- ifelse(data$pcor > 0, FALSE, TRUE)
   			 		edges$dashes <- dashes
    				edges$values <- data$pcor
    				edges$title <- paste(edges$node1, "-", 
    								edges$node2, " : ", edges$values, sep="")
   			 	} else if(calc_type %in% "multibipartite_ggm") {
   			 		dashes <- ifelse(data$coeffs.1 > 0 & data$coeffs.2 > 0, FALSE, TRUE)
   			 		edges$dashes <- dashes
    				edges$values <- (data$coeffs.1 + data$coeffs.2)/2 
    				edges$title <- paste("coeff1: ", 
    					data$coeffs.1, "<br /> coeff2: ", data$coeffs.2, sep=" ")
   			 	} else if(calc_type %in% "temporal_network") {
   			 		dashes <- ifelse(data$coeffs > 0, FALSE, TRUE)
   			 		edges$values <- data$coeffs
   			 		edges$dashes <- dashes
   			 		edges$title <- paste(edges$node1, "-", edges$node2, " : ", edges$values, sep="")
   			 		edges$arrows <- rep("from", each=length(edges$dashes))
   			 	}
   			 	if(length(grep("calc_|mod_merge_results|add_result|meta_", object@results[[length(object@results)]]$functions_applied))==1) {
					object@results[[length(object@results)+1]] <- list(functions_applied=c(), 
						plot_data=list(network=list(node=nodes, edge=edges)),
						information=list(calc_type=calc_type, calc_info=calc_info), plots=list())
				} else {
					object@results[[length(object@results)]]$plot_data$network$node <- nodes
					object@results[[length(object@results)]]$plot_data$network$edge <- edges
					object@results[[length(object@results)]]$information$calc_type <- calc_type
					object@results[[length(object@results)]]$information$calc_info <- calc_info
					object@results[[length(object@results)]]$plots <- list()
				}
			} else if(length(grep("merged", calc_type)==1)) {
				if(length(grep("calc_|mod_merge_results|add_result|meta_", object@results[[length(object@results)]]$functions_applied))==1) {
					object@results[[length(object@results)+1]] <- list(functions_applied=c(), 
						plot_data=list(merged_network=data),
						information=list(calc_type=calc_type, calc_info=calc_info), plots=list())
				} else {
					object@results[[length(object@results)]]$plot_data[["merged_network"]] <- data
					object@results[[length(object@results)]]$information$calc_type <- calc_type
					object@results[[length(object@results)]]$information$calc_info <- calc_info
					object@results[[length(object@results)]]$plots <- list()
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
						dummy_metadata$id <- NULL
						return(cbind.data.frame(data[[x]], dummy_metadata))
					})
				}
				if(length(grep("calc_|mod_merge_results|add_result|meta_", object@results[[length(object@results)]]$functions_applied))==1) {
					object@results[[length(object@results)+1]] <- list(functions_applied=c(), plot_data=plot_data,
							information=list(calc_type=calc_type, calc_info=calc_info), plots=list())
				} else {
					object@results[[length(object@results)]]$plot_data <- plot_data
					object@results[[length(object@results)]]$information$calc_type <- calc_type
					object@results[[length(object@results)]]$information$calc_info <- calc_info
					object@results[[length(object@results)]]$plots <- list()
				}
			}
			names(object@results)[length(object@results)] <- name
			return(object)
})


