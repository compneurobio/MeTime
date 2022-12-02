#' Setting up standard wrapper for network plots from visNetwork for a metime_plotter object. 
#' @description plot function for metime_plotter object with different inputs to specialize plots. Used for all calc outputs.
#' @param object S4 object of class metime_plotter
#' @param title character/string that is the title of the graph output
#' @param layout_by character to define the layout style to be used
#' @export
setGeneric("viz_plotter_visNetwork", function(object, title, layout_by) standardGeneric("viz_plotter_visNetwork"))
setMethod("viz_plotter_visNetwork", "metime_plotter", function(object, title, layout_by) { 
		metadata <- object@plot_data[["metadata"]]
		node_list <- object@plot_data[["node"]]
		edge_list <- object@plot_data[["edge"]]
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
            graph <- visNetwork::visNetwork(nodes=object@plot_data[["node"]], edges=object@plot_data[["edge"]], main=title) %>%
                    visNetwork::visIgraphLayout(layout=layout_by, physics = F, smooth = F) %>%
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
        	graph <- visNetwork::visNetwork(nodes=object@plot_data[["node"]], edges=object@plot_data[["edge"]], main=title) %>%
                    visNetwork::visIgraphLayout(layout=layout_by, physics = F, smooth = F) %>%
                    visNetwork::visPhysics(stabilization = FALSE) %>%
                    visNetwork::visLegend(useGroups = T) %>% 
                    visNetwork::visNodes(borderWidth = 3) %>%
                    visNetwork::visEdges(smooth = FALSE, shadow = TRUE) %>%
                    visNetwork::visOptions(highlightNearest = list(enabled=T, hover=T), nodesIdSelection = T, selectedBy = "group") %>%
                    visNetwork::visInteraction(navigationButtons = T) %>%
                    visNetwork::visExport(type="pdf", name=title, float="right")		
        }
        object@plot[[1]] <- graph
        out <- object
        return(out)
	})





