#' Function to extract plot data into a csv
#' @description extracts information from plotter object and saves it as a csv
#' @examples see examples here
#' Network : save_plotter_object(object, out=c("edge.csv", "node.csv", "meta.csv"))
#' Others : save_plotter_object(object, out="outfile.csv")
#' @param object An object of class metime_plotter
#' @param out Character to specify path of the output file or character vector in case of visNetwork
#' @return saves the data into a csv and returns nothing
#' @export
setGeneric("save_plotter_object", function(object, out) standardGeneric("save_plotter_object"))
setMethod("save_plotter_object", "metime_plotter", function(object, out) {
			if(object@style %in% "visNetwork") {
				node <- plot_data[["node"]] 
				edge <-	plot_data[["edge"]] 
				meta <-	plot_data[["metadata"]]
				write.csv(node, out[2])
				write.csv(edge, out[1])
				write.csv(meta, out[3]) 
			} else {
				plot_data <- object@plot_data
				calc_type <- rep(object@calc_type, each=length(plot_data[,1]))
				calc_info <- rep(object@calc_info, each=length(plot_data[,1]))
				plot_data <- cbind(plot_data, calc_info, calc_type)
				write.csv(plot_data, out)
			}
	}) 


#' Function to save interactive plots
#' @description extracts plot from plotter object and saves it as a widget
#' @examples save_plot_from_plotter(object)
#' @param object An object of class metime_plotter
#' @param out Character to specify path of the output file to save the widget in
#' @return saves the plot and returns nothing
#' @export
setGeneric("save_plot_from_plotter", function(object, out) standardGeneric("save_plot_from_plotter"))
setMethod("save_plot_from_plotter", "metime_plotter", function(object, out) {
			stopifnot(length(object@plot) == length(out))
			plots <- object@plot
			for(i in 1:length(out)) {
				saveWidget(plots[[i]], out[i])
			}
	})

#' Function to extract analyser object data into a csv
#' @description extracts information from analyser object and saves it as a csv
#' @examples see examples here
#' save_analyser_object(object, which_data="dataset")
#' @param object An object of class metime_plotter
#' @param which_data Character to specify the dataset
#' @return saves the data in the working directory as a csv and returns nothing
#' @export
setGeneric("save_analyser_object", function(object, which_data) standardGeneric("save_analyser_object"))
setMethod("save_analyser_object", "metime_analyser", function(object, which_data) {
			list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
			list_of_col_data <- object@list_of_col_data[names(object@list_of_col_data) %in% which_data]
			list_of_row_data <- object@list_of_row_data[names(object@list_of_row_data) %in% which_data]
			for(i in 1:length(which_data)) {
				write.csv(list_of_data[[i]], paste(which_data[i], ".csv", sep=""))
				write.csv(list_of_col_data[[i]], paste(which_data[i], "_col.csv", sep=""))
				write.csv(list_of_row_data[[i]], paste(which_data[i], "_row.csv", sep=""))		
			}
	}) 
