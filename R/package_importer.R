usethis::use_import_from("rstatix", fun="t_test")
usethis::use_package("tidyverse", type="depends")
#usethis::use_package("igraph", type="Imports")
#usethis::use_package("plyr", type="Imports")
usethis::use_package("DT", type="Imports")
usethis::use_import_from("M3C", fun="tsne")
usethis::use_import_from("umap", fun="umap")
usethis::use_package("plotly", type="Imports")
usethis::use_package("visNetwork", type="Imports")
usethis::use_package("RColorBrewer", type="Imports")
usethis::use_package("GeneNet", type="Imports")
usethis::use_package("longitudinal", type="Imports")
usethis::use_package("glmnet", type="Imports")
usethis::use_package("Boruta", type="Imports")
usethis::use_import_from("nnet", fun="class.ind") # get that manually 
usethis::use_import_from("abind", fun="abind") # get that manually
usethis::use_import_from("multiway", fun="parafac") # get that manually
#usethis::use_import_from("DT", fun="datatable")
usethis::use_import_from("xlsx", fun="write.xlsx") # get that manually, or check for base function or in tidyverse
usethis::use_package("visNetwork", type="Imports")
usethis::use_package("mgcv", type="Imports")
usethis::use_import_from("lmerTest", fun="lmer")
usethis::use_import_from("parallel", fun="mclapply")
#usethis::use_import_from("reshape2", fun="melt")
#usethis::use_import_from("data.table", fun="setDT") 
usethis::use_package("MatrixEQTL", type="Imports")
usethis::use_package("WGCNA", type="Imports") # change it to functions
usethis::use_import_from("dynamicTreeCut", fun="cutreeDynamic")


#' Validity function to check if the object is valid or not
#' @description Function to check the validity of metime_analyser 
#' @param object An S4 object of class metime_analyser
#' @examples validObject(object)
#' @returns logical suggesting if the object is intact or not
#' @export
validity <- function(object) {
		out <- TRUE
		if(!check_rownames_and_colnames(object)) out <- FALSE 
		if(!check_ids_and_classes(object)) out <- FALSE
		if(!check_results(object)) out <- FALSE
		return(out) 
	}

#creating reference metime-analyser class that creates an object with full data

#' Constructor to generate an object of class metime_analyser. 
#' contains slots - list_of_data: For the list of all data matrices.
#'				  - list_of_col_data: list of all the col data files in the same order.
#' 				  - list_of_row_data: list of all the row data files in the same order.
#' 				  - annotations: list with phenotype and medication. Each of which is character that represents 
#'									the name of the aforementioned dataset types.  	
#' 	
#' @rdname metime_analyser
#' @export 
setClass("metime_analyser", slots=list(list_of_data="list", list_of_col_data="list", list_of_row_data="list", 
								 annotations="list", results="list"), validity=validity) 



# add a function to generate Rmarkdown directly from the analyser object - viz_analyser()
# create another function which add the function applied on the analyser_object - ongoing
# Try to improve the code of all functions and then Matthias will check it

#' Setting a plotting method for the metime_analyser class
#' @description Function to plot results of a certain calculation 
#' @param object An S4 object of class metime_analyser
#' @param results_index Index/name of the results to be plotted
#' @param interactive logical. Set TRUE for interactive plot
#' @param ... other parameters to pass color, fill, strat, viz(character vector with colnames for interactive)
#'  
#' @return plots for a certain set of results
#' @export
setMethod("plot", "metime_analyser", function(x, results_index, interactive, ...) {
			add <- list(...)
			results <- x@results[[results_index]]
			if(!all(names(plot_data) %in% c("node", "edge", "metadata"))) {
				plots <- lapply(seq_along(results$plot_data), function(x) {
					if(results$information$calc_type[x] %in% "CI_metabolite" |
						results$information$calc_type[x] %in% "CI_metabotype") {
						plot <- ggplot(results$plot_data[[x]], aes(x=x, y=ci)) + 
								geom_point(aes_string(color=add$color, 
								shape=add$shape)) + facet_wrap(add$strats) + theme_classic()
						if(interactive) {
							plot <- plotly::ggplotly(plot, width=800, height=800)
							for(i in seq_along(plot$x$data)) {
								x <- plot$x$data[[i]]$x
								data <- results$plot_data[[x]]
								data <- data[data[ ,data$x] %in% x, ] 
								plot$x$data[[j]]$text <- get_text_for_plot(data=data, 
										colnames=add$viz)
							}
							return(plot)
						} else {
							return(plot)
						}
					} else if(results$information$calc_type[x] %in% "pairwise_distance" |
						results$information$calc_type[x] %in% "pairwise_correlation") {
						plot <- ggplot(results$plot_data[[x]], aes(x=row,y=column)) + 
							geom_tile(aes_string(fill=dist)) +
							facet_wrap(add$strats) + theme_classic()
						if(interactive) {
							plot <- plotly::ggplotly(plot, width=800, height=800)
							for(i in seq_along(plot$x$data)) {
								x <- plot$x$data[[i]]$x
								data <- results$plot_data[[x]]
								data <- data[data[ ,data$row] %in% x, ] 
								plot$x$data[[j]]$text <- get_text_for_plot(data=data, 
									colnames=add$viz)
							}
							return(plot)
						} else {
							return(plot)
						}
					} else if(results$information$calc_type[x] %in% "PCA") {
						plot <- ggplot(results$plot_data[[x]], aes(x=PC1, y=PC2)) + 
								geom_point(aes_string(color=add$color, 
								shape=add$shape)) + facet_wrap(add$strats) + theme_classic()
						if(interactive) {
							plot <- plotly::ggplotly(plot, width=800, height=800)
							for(i in seq_along(plot$x$data)) {
								x <- plot$x$data[[i]]$x
								data <- results$plot_data[[x]]
								data <- data[data[ ,data$PC1] %in% x, ] 
								plot$x$data[[j]]$text <- get_text_for_plot(data=data, colnames=add$viz)
							}
							return(plot)
						} else {
							return(plot)
						}
					} else if(results$information$calc_type[x] %in% "UMAP") {
						plot <- ggplot(results$plot_data[[x]], aes(x=UMAP1, y=UMAP2)) + 
								geom_point(aes_string(color=add$color, 
								shape=add$shape)) + facet_wrap(add$strats) + theme_classic()
						if(interactive) {
							plot <- plotly::ggplotly(plot, width=800, height=800)
							for(i in seq_along(plot$x$data)) {
								x <- plot$x$data[[i]]$x
								data <- results$plot_data[[x]]
								data <- data[data[ ,data$UMAP1] %in% x, ] 
								plot$x$data[[j]]$text <- get_text_for_plot(data=data, colnames=add$viz)
							}
							return(plot)
						} else {
							return(plot)
						}
					} else if(results$information$calc_type[x] %in% "tSNE") {
						plot <- ggplot(results$plot_data[[x]], aes(x=X1, y=X2)) + 
								geom_point(aes_string(color=add$color, 
								shape=add$shape)) + facet_wrap(add$strats) + theme_classic()
						if(interactive) {
							plot <- plotly::ggplotly(plot, width=800, height=800)
							for(i in seq_along(plot$x$data)) {
								x <- plot$x$data[[i]]$x
								data <- results$plot_data[[x]]
								data <- data[data[ ,data$X1] %in% x, ] 
								plot$x$data[[j]]$text <- get_text_for_plot(data=data, colnames=add$viz)
							}
							return(plot)
						} else {
							return(plot)
						}
					}
					
				})
				return(plots)
			} else {
				metadata <- results$plot_data$metadata
				node_list <- results$plot_data$node
				edge_list <- results$plot_data$edge
			#Choosing colors and shapes for visualization
        		if(is.null(node_list$color)) {
            		shapes <- c("square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond")
            		colors_for_nodes <- node_list$group
            		colors_code <- as.data.frame(cbind(unique(node_list$group), get_palette(length(unique(node_list$group)))))
            		colnames(colors_code) <- c("group", "color")
            		colors_new <- c()
            		for(i in 1:length(colors_for_nodes)) {
            			colors_new[i] <- colors_code$color[colors_code$group %in% colors_for_nodes[i]]
            		}
            		color <- colors_new
            		node_list <- as.data.frame(cbind(node_list, color))
        		} else {
            		shapes <- c("square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond")
        		}
        		if(length(unique(metadata$class)) > 1) {
        			classes <- unique(metadata$class)
            		groups <-  unique(metadata$group)
            		graph <- visNetwork::visNetwork(nodes=node_list, edges=edge_list) %>%
                    	visNetwork::visIgraphLayout(layout="layout.fruchterman.reingold", physics = F, smooth = F) %>%
                    	visNetwork::visPhysics(stabilization = FALSE)
            		for(i in 1:length(classes)) {
                		graph <- graph %>% visNetwork::visGroups(graph=graph, groupname = classes[i], shape = shapes[i])
            		}
            		graph <- graph %>% visNetwork::visLegend(useGroups = T) %>% 
                    	visNetwork::visNodes(borderWidth = 3, color=list(background=colors_for_nodes)) %>%
                    	visNetwork::visEdges(smooth = FALSE, shadow = TRUE) %>%
                        visNetwork::visOptions(highlightNearest = list(enabled=T, hover=T), nodesIdSelection = T, selectedBy = "group") %>%
                        visNetwork::visInteraction(navigationButtons = T) %>%
                       	visNetwork::visExport(type="pdf", name=title, float="right")
        		} else {
        			graph <- visNetwork::visNetwork(nodes=node_list, edges=edge_title) %>%
                    	visNetwork::visIgraphLayout(layout="layout.fruchterman.reingold", physics = F, smooth = F) %>%
                    	visNetwork::visPhysics(stabilization = FALSE) %>%
                    	visNetwork::visLegend(useGroups = T) %>% 
                    	visNetwork::visNodes(borderWidth = 3) %>%
                    	visNetwork::visEdges(smooth = FALSE, shadow = TRUE) %>%
                    	visNetwork::visOptions(highlightNearest = list(enabled=T, hover=T), nodesIdSelection = T, selectedBy = "group") %>%
                    	visNetwork::visInteraction(navigationButtons = T) %>%
                    	visNetwork::visExport(type="pdf", name=title, float="right")		
        		}
        		return(graph)
        	}
	})

#' Setting new print definition for the metime_analyser object
#' @description function to see the structure of metime_analyser object
#' @param object S4 object of class metime_analyser
#' @examples structure(object)
#' @return structure of the S4 object
#' @export
setMethod("show", "metime_analyser", function(object) {
		out <- lapply(seq_along(object@list_of_data), function(i) {
			if(names(object@list_of_data)[i] %in% object@annotations$phenotype) {
				cat(paste("Dataset:", names(object@list_of_data)[i], "with", 
				dim(object@list_of_data[[i]])[1], "samples", "and", dim(object@list_of_data[[i]])[2], "phenotypes",sep=" "))
				cat("\n")
			} else if(names(object@list_of_data)[i] %in% object@annotations$medication) {
				cat(paste("Dataset:", names(object@list_of_data)[i], "with", 
				dim(object@list_of_data[[i]])[1], "samples", "and", dim(object@list_of_data[[i]])[2], "medicines", sep=" "))
				cat("\n")
			} else {
				cat(paste("Dataset:", names(object@list_of_data)[i], "with", 
				dim(object@list_of_data[[i]])[1], "samples", "and", dim(object@list_of_data[[i]])[2], "metabolites",sep=" "))
				cat("\n")
			}
			return(NULL)
		})
	})
