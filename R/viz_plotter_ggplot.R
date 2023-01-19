
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
					object@plot[[i]] <- object@plot[[i]] + geom_point(aes_string(color=aesthetics[[i]]$color, 
										shape=aesthetics[[i]]$shape)) + 
										facet_wrap(aesthetics[[i]]$strats) +
										theme_classic()
				} else if(object@plot_type[i] %in% "heatmap") {
					object@plot[[i]] <- object@plot[[i]] + geom_tile(aes_string(fill=aesthetics[[i]]$fill)) + 
										facet_wrap(aesthetics[[i]]$strats) +
										theme_classic()
				} else if(object@plot_type[i] %in% "line") {
					object@plot[[i]] <- object@plot[[i]] + geom_line(aes_string(color=aesthetics[[i]]$color)) + 
										facet_wrap(aesthetics[[i]]$strats) + 
										geom_point(aes_string(color=aesthetics[[i]]$color)) +
										theme_classic()
				} else if(object@plot_type[i] %in% "box") {
					object@plot[[i]] <- object@plot[[i]] + geom_boxplot(aes_string(color=aesthetics[[i]]$color), 
										outlier.colour="red", outlier.shape=8, outlier.size=4) + 
										facet_wrap(aesthetics[[i]]$strats) +
										theme_classic() + stat_summary(fun.y=mean, geom="point", shape=23, size=4)
				} else if(object@plot_type[i] %in% "forest") {
					object@plot[[i]] <- object@plot[[i]] + geom_pointrange(aes_string(ymin=aesthetics[[i]]$lower, 
										ymax=aesthetics[[i]]$upper, color=aesthetics[[i]]$color)) + 
										coord_flip() + facet_wrap(aesthetics[[i]]$strats) +
										theme_classic() 
				} else if(object@plot_type[i] %in% "QQ") {
					object@plot[[i]] <- object@plot[[i]] + stat_qq(aes_string(sample=aesthetics[[i]]$sample, 
															color=aesthetics[[i]]$color, 
															shape=aesthetics[[i]]$shape)) + 
							facet_wrap(aesthetics[[i]]$strats) + theme_classic()	
				} else if(object@plot_type[i] %in% "bar") {
					object@plot[[i]] <- object@plot[[i]] + geom_bar(aes_string(fill=aesthetics[[i]]$fill), 
										stat="identity", position=aesthetics[[i]]$position) + 
										facet_wrap(aesthetics[[i]]$strats) + theme_classic()
				} else {
					stop("input for plot_type is either wrong or such a plot is not available")
				}

				if(interactive) {
					object@plot[[i]] <- plotly::ggplotly(object@plot[[i]], width = 800, height = 600)
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

