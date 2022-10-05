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
								 annotations="list")) 


#' creating metime_plotter class that converts calculations and metadata as a plotable object to parse 
#' into viz_plotter
#' Contains slots - plot_data: Dataframe with plotting data and metadata for visualization
#' 				  - plot: ggplot(), circos() or visNetwork() object with predefined aesthetics 
#'                - calc_type: A vector to specify type of calculation - will be used for comp_ functions
#'                - calc_info: string to define the information about calculation
#' 				  - plot_type: A character vector to define the type of plots that are needed.
#' @rdname metime_plotter
#' @export
setClass("metime_plotter", slots=list(plot_data="list", plot="list", calc_type="character", calc_info="character", plot_type="character", style="character"))


#' Function for Plotting distributions of phenotypic variables 
#' @description A method to be applied onto s4 object so as to obtain distributions of various phenotypic variables
#' @examples # extracting distribuiton of Age from dataset1
#' plot <- viz_distribution_plotter(object, colname="Age", which_data="dataset1", strats="additional columns for facet wrapping")
#' @param object An object of class metime_analyser
#' @param colname Name of the variable whose distribution is of interest
#' @param which_data Name of the dataset from which the samples will be extracted
#' @param strats Character vector with colnames that are to be used for stratification
#' @param phenotype Logical. If true data will be collected from phenotype_data matrix else from row data 
#' @return a list with either 1) density plot, mean table acc to timepoint and variable type or 
#'								2) bar plot, line plot, and variable type
#' @export
setGeneric("viz_distribution_plotter", function(object, colname, which_data, strats, phenotype) standardGeneric("viz_distribution_plotter") )

setMethod("viz_distribution_plotter", "metime_analyser",function(object, colname, which_data, strats=NULL, phenotype) {
	data <- object@list_of_data[names(object@list_of_data) %in% which_data]
	data <- data[[1]]
	if(phenotype) {
		phenotype <- object@list_of_data[[object@annotations$phenotype]]
	} else {
		phenotype <- object@list_of_row_data[[which_data]]
	}
	var_type <- c()
	vec <- phenotype[,colname]
	vec <- na.omit(vec)
	if(is.character(vec)==TRUE) var_type <- "bar"
	if(is.numeric(vec)==TRUE) {
		if(length(unique(vec)) <= 12) {
				var_type <- "bar"
		} else {
				var_type <- "density"
		}
	}
	if(var_type %in% "bar") {
		phenotype <- phenotype[rownames(phenotype) %in% rownames(data),]
		vec <- phenotype[, colname]
        palette_line <- get_palette(length(unique(vec)))
		palette_timepoints <- get_palette(30)[c(16,2,3,6,7,8,10,13,19,20,23,27,30)]
		timepoints <- unlist(lapply(strsplit(rownames(phenotype), split="_"), function(x) return(x[2])))
		levels <- sort(unique(as.numeric(unlist(lapply(strsplit(timepoints, split="t"), function(x) return(x[2]))))))
		timepoints <- factor(timepoints, levels=paste("t",levels,sep=""))
		plot_data <- data.frame(vec=as.character(vec), timepoints=timepoints)
		plot_data <- table(plot_data)
		plot_data <- reshape2::melt(plot_data)
		colnames(plot_data) <- c(colname, "Timepoints", "Frequency")
		plot_data[,colname] <- as.factor(plot_data[,colname])
		bar_plot <- ggplot(data=plot_data, aes_string(x=colname, y="Frequency", fill="Timepoints")) +
					geom_bar(stat="identity") + scale_fill_manual(values=palette_timepoints)+ theme_classic()
		line_plot <- ggplot(plot_data, aes_string(x="Timepoints", y="Frequency", group=colname)) + 
					 geom_line(aes_string(color=colname)) + geom_point(aes_string(color=colname)) + scale_color_manual(values=palette_line) + theme_classic()
		return(out=list(bar_plot=bar_plot, line_plot=line_plot, var_type=var_type))

	} else {
		phenotype <- phenotype[rownames(phenotype) %in% rownames(data),]
		vec <- phenotype[,colname]
		palette_timepoints <- get_palette(30)[c(16,2,3,6,7,8,10,13,19,20,23,27,30)]
		timepoints <- unlist(lapply(strsplit(rownames(phenotype), split="_"), function(x) return(x[2])))
		if(is.null(strats)) {
			plot_data <- as.data.frame(cbind(as.character(vec), timepoints))
			plot_data[,1] <- as.numeric(plot_data[,1])
			levels <- sort(unique(as.numeric(unlist(lapply(strsplit(timepoints, split="t"), function(x) return(x[2]))))))
			plot_data$timepoints <- factor(timepoints, levels=paste("t", levels, sep=""))
			colnames(plot_data) <- c(colname, "Timepoints")
			plot_data <- na.omit(plot_data)
			mu <- as.data.frame(aggregate(plot_data[,1], list(Timepoints=plot_data[,2]), FUN=mean))
			mu$Timepoints <- factor(mu$Timepoints, levels=paste("t", levels, sep=""))
			colnames(mu)[2] <- "mean"
			density_plot <- ggplot(plot_data, aes_string(x=colname, fill="Timepoints")) + geom_density(alpha=0.3) + 
								geom_vline(data=mu, aes(xintercept=mean, color=Timepoints), linetype="dashed") +
			 					scale_fill_manual(values=palette_timepoints) + theme_classic()
		} else {
			strats_df <- phenotype[,strats]
			plot_data <- as.data.frame(cbind(as.character(vec), timepoints, strats_df))
			plot_data[,1] <- as.numeric(plot_data[,1])
			levels <- sort(unique(as.numeric(unlist(lapply(strsplit(timepoints, split="t"), function(x) return(x[2]))))))
			plot_data$timepoints <- factor(timepoints, levels=paste("t", levels, sep=""))
			colnames(plot_data) <- c(colname, "Timepoints", strats)
			plot_data <- na.omit(plot_data)
			mu <- as.data.frame(aggregate(plot_data[,1], list(Timepoints=plot_data[,2]), FUN=mean))
			mu$Timepoints <- factor(mu$Timepoints, levels=paste("t", levels, sep=""))
			colnames(mu)[2] <- "mean"
			density_plot <- ggplot(plot_data, aes_string(x=colname, fill="Timepoints")) + geom_density(alpha=0.3) + 
								geom_vline(data=mu, aes(xintercept=mean, color=Timepoints), linetype="dashed") +
			 					scale_fill_manual(values=palette_timepoints) + theme_classic() + facet_wrap(strats)
		}
		
		#render table for mean values
		out <- list(table=mu, density_plot=density_plot, var_type=var_type)
		return(out)
	}
})

# #' Function to dot plot any kind of dot_plotter including for dimensionality reduction
# #' @description General function to be implemented on data_list that is obtained after applying a dimensionality reduction method
# #' @param data_list list obtained after applying calc_dimensionality_reduction() on metime_analyse object
# #' @param metadata_list list obtained after applying get_metadata_for_plotting() on metime_analyse object
# #' @param axes_labels character vector to specify the labels of the axes in the order x and y.
# #' @param title_metabs character to specify the title of the plot of metabolites
# #' @param title_metabs character to specify the title of the plot of metabolites
# #' @return a list with both the plots of samples and metabolites. Can be accessed by using ".$samples" and ".$metabs"
# #' @export
# viz_dimensionality_reduction <- function(data_list, metadata_list, axes_labels, title_metabs, title_samples) {
# 	palette <- get_palette(30)
# 	if(is.null(metadata_list)) {
# 		return(list(metabs=plot(data_list$metabs), samples=plot(data_list$samples)))
# 	} else {
# 		plot <- list()
# 		data_list <- lapply(data_list, function(x) return(x[sort(rownames(x)),]))
# 		plot_data_metabs <- as.data.frame(cbind(data_list$metabs, metadata_list$metabs))
# 		plot_data_samples <- as.data.frame(cbind(data_list$samples, metadata_list$samples))
# 		plot_metabs <- ggplot(plot_data_metabs, aes_string(x=colnames(plot_data_metabs)[1], y=colnames(plot_data_metabs)[2], 
# 							color="groups", shape="class", text=paste("metab_name :", plot_data_metabs$name, sep="")))  
# 						+ geom_point() + labs(x=axes_labels[1], y=axes_labels[2], subtitle=title_metabs) 
# 		plot_metabs <- ggplotly(metabs)
# 		timepoints <- unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[2])))
# 		levels <- unique(sort(as.numeric(unlist(lapply(strsplit(timepoints, split="t", fixed=TRUE), function(x) return(x[2]))))))
# 		plot_data_samples$timepoints <- factor(timepoints, levels=paste("t", levels, sep=""))
# 		plot_samples <- ggplot(plot_data_samples, aes_string(x=colnames(plot_data_samples)[1], y=colnames(plot_data_samples)[2], 
# 						color="timepoints", 
# 						text=get_text_for_plot(data=plot_data_samples, colnames=colnames(plot_data_samples)[3:length(colnames(plot_data_samples))]))) +
# 						geom_point() + labs(x=axes_labels[1], y=axes_labels[2], subtitle=title_samples)
# 		plot_samples <- ggplotly(plot_samples)
# 		return(list(metabs=plot_metabs, samples=plot_samples))
# 	}
# }
 
#viz_conversation_index()

#Create an S4 class for plotting
#analysis - needs ids or names, col_id data, row_data, col_data

#set slots in parameter_list and keep it as an inheritance

#make_plot_object function to be written

#' Setting up standard wrapper for all ggplot plots for a metime_plotter object. 
#' @description plot function for metime_plotter object with different inputs to specialize plots. Used for all calc outputs.
#' @param object S4 object of class metime_plotter
#' @param aesthetics list for aesthetics. eg: list(list(x="colname",y="colname",color="colname", shape="colname"), list(...)) for "dot" plot and "heatmap"
#' plot, for heatmap: list(x="colname", y="colname", fill="colname"). Additionally two other character vectors are allowed namely .$viz and .$strats for text
#' and for facet wrapping. 
#' @param interactive Flag option(Logical). If set to TRUE will generate an interactive plot else will generate a normal ggplot
#' @return metime_plotter object with updated plot
#' @export
setGeneric("viz_plotter_ggplot", function(object, aesthetics, interactive) standardGeneric("viz_plotter_ggplot"))
setMethod("viz_plotter_ggplot", "metime_plotter", function(object, aesthetics, interactive) {
			stopifnot(object@plot_type %in% c("dot", "heatmap", "line", "box", "bar", "forest", "QQ"))
			for(i in 1:length(object@plot_type)) {
				if(object@plot_type[i] %in% "dot") {
					object@plot[[i]] <- object@plot[[i]] + geom_point(aes_string(x=aesthetics[[i]]$x, y=aesthetics[[i]]$y, 
										color=aesthetics[[i]]$color, shape=aesthetics[[i]]$shape)) + 
										facet_wrap(aesthetics[[i]]$strats) +
										theme_classic()
				} else if(object@plot_type[i] %in% "heatmap") {
					object@plot[[i]] <- object@plot[[i]] + geom_tile(aes_string(x=aesthetics[[i]]$x, y=aesthetics[[i]]$y, 
										fill=aesthetics[[i]]$fill)) + facet_wrap(aesthetics[[i]]$strats) +
										theme_classic()
				} else if(object@plot_type[i] %in% "line") {
					object@plot[[i]] <- object@plot[[i]] + geom_line(aes_string(x=aesthetics[[i]]$x, y=aesthetics[[i]]$y, 
										color=aesthetics[[i]]$color)) + facet_wrap(aesthetics[[i]]$strats) + geom_point(aes_string(color=aesthetics[[i]]$color)) +
										theme_classic()
				} else if(object@plot_type[i] %in% "box") {
					object@plot[[i]] <- object@plot[[i]] + geom_boxplot(aes_string(x=aesthetics[[i]]$x, y=aesthetics[[i]]$y, 
										color=aesthetics[[i]]$color), outlier.colour="red", outlier.shape=8, outlier.size=4) + facet_wrap(aesthetics[[i]]$strats) +
										theme_classic() + stat_summary(fun.y=mean, geom="point", shape=23, size=4)
				} else if(object@plot_type[i] %in% "forest") {
					object@plot[[i]] <- object@plot[[i]] + geom_pointrange(aes_string(x=aesthetics[[i]]$label, 
																	y=aesthetics[[i]]$mean, ymin=aesthetics[[i]]$lower, 
																	ymax=aesthetics[[i]]$upper, color=aesthetics[[i]]$color)) + 
										coord_flip() + facet_wrap(aesthetics[[i]]$strats) +
										theme_classic() 
				} else if(object@plot_type[i] %in% "QQ") {
					object@plot[[i]] <- object@plot[[i]] + stat_qq(aes_string(sample=aesthetics[[i]]$sample, color=aesthetics[[i]]$color, 
															shape=aesthetics[[i]]$shape)) + facet_wrap(aesthetics[[i]]$strats) + theme_classic()	
				} else if(object@plot_type[i] %in% "bar") {
					object@plot[[i]] <- object@plot[[i]] + geom_bar(aes_string(x=aesthetics[[i]]$x, y=aesthetics[[i]]$y, 
										fill=aesthetics[[i]]$fill), stat="identity", position=aesthetics[[i]]$postion) + 
										facet_wrap(aesthetics[[i]]$strats) + theme_classic()
				} else {
					stop("input for plot_type is either wrong or such a plot is not available")
				}

				if(interactive) {
					object@plot[[i]] <- ggplotly(object@plot[[i]], width = 800, height = 600)
					for(j in 1:length(object@plot[[i]]$x$data)) {
						x <- object@plot[[i]]$x$data[[j]]$x
						data <- object@plot_data[[1]]
						data <- data[data[ ,aesthetics[[i]]$x] %in% x, ] 
						object@plot[[i]]$x$data[[j]]$text <- get_text_for_plot(data=data, colnames=aesthetics[[i]]$viz)
					}
				}
			}
			out <- object
			return(out)
			
	})

#lol <- viz_plotter_ggplot(object=plotter_object, aesthetics=list(list(x="PC1",y="PC2",color="group", shape="class"))) 

#' Setting up standard wrapper for all circos plots for a metime_plotter object. 
#' @description plot function for metime_plotter object with different inputs to specialize plots. Used for all calc outputs.
#' @param object S4 object of class metime_plotter
#' @param aesthetics list for aesthetics. eg: list(list(x="colname",y="colname",color="colname", shape="colname"), list(...)) for "dot" plot and "heatmap"
#' plot, for heatmap: list(x="colname", y="colname", fill="colname"). Additionally two other character vectors are allowed namely .$vis and .$strats for text
#' and for facet wrapping.
#' @export
setGeneric("viz_plotter_circos", function(object, aesthetics, outfile, layout_by) standardGeneric("viz_plotter_circos")) 
setMethod("viz_plotter_circos", "metime_plotter", function(object, aesthetics, outfile, layout_by) {
			pdffile <- pdf(outfile)

			dev.off()
	}) 

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
        if(is.null(node_list$colors)) {
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
            graph <- visNetwork(nodes=object@plot_data[["node"]], edges=object@plot_data[["edge"]], main=title) %>%
                    visIgraphLayout(layout=layout_by, physics = F, smooth = F) %>%
                    visPhysics(stabilization = FALSE)
            for(i in 1:length(classes)) {
                graph <- graph %>% visGroups(graph=graph, groupname = classes[i], shape = shapes[i])
            }
            graph <- graph %>% visLegend(useGroups = T) %>% 
                    visNodes(borderWidth = 3, color=list(background=colors_for_nodes)) %>%
                    visEdges(smooth = FALSE, shadow = TRUE) %>%
                        visOptions(highlightNearest = list(enabled=T, hover=T), nodesIdSelection = T, selectedBy = "group") %>%
                        visInteraction(navigationButtons = T) %>%
                       	visExport(type="pdf", name=title, float="right")%>%
 						visConfigure(enabled = TRUE)
        } else {
        	graph <- visNetwork(nodes=object@plot_data[["node"]], edges=object@plot_data[["edge"]], main=title) %>%
                    visIgraphLayout(layout=layout_by, physics = F, smooth = F) %>%
                    visPhysics(stabilization = FALSE) %>%
                    visLegend(useGroups = T) %>% 
                    visNodes(borderWidth = 3) %>%
                    visEdges(smooth = FALSE, shadow = TRUE) %>%
                    visOptions(highlightNearest = list(enabled=T, hover=T), nodesIdSelection = T, selectedBy = "group") %>%
                    visInteraction(navigationButtons = T) %>%
                    visExport(type="pdf", name=title, float="right")%>%
  					visConfigure(enabled = TRUE)
        }
        object@plot[[1]] <- graph
        out <- object
        return(out)
	})







#Update layout by in the function -done



