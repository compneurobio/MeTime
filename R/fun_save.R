#' Function to extract plot data into a csv
#' @description extracts information from plotter object and saves it as a csv
#' @examples see examples here and maintain the order
#' If type is xlsx one out file is enough and the network data will be stored in different sheets
#' Network : save_plotter_object(object, out=c("edge", "node", "meta"), type="csv")
#' Others : save_plotter_object(object, out="outfile", type="tsv")
#' @param object An object of class metime_plotter
#' @param out Character to specify path of the output file or character vector in case of visNetwork
#' @param type character to define outfile type that is "csv", "xlsx" or "tsv"
#' @return saves the data into a csv and returns nothing
#' @export
setGeneric("save_plotter_object", function(object, out, type) standardGeneric("save_plotter_object"))
setMethod("save_plotter_object", "metime_plotter", function(object, out, type) {
			stopifnot(type %in% c("csv", "xlsx", "tsv"))
			if(object@style %in% "visNetwork") {
				node <- plot_data[["node"]] 
				edge <-	plot_data[["edge"]] 
				meta <-	plot_data[["metadata"]]
				if(type %in% "csv") {
					write.csv(node, paste(out[2], ".csv", sep=""), row.names=FALSE)
					write.csv(edge, paste(out[1], ".csv", sep=""), row.names=FALSE)
					write.csv(meta, paste(out[3], ".csv", sep=""), row.names=FALSE)
				} else if(type %in% "tsv") {
					write.table(node, file=paste(out[2], ".tsv", sep=""), quote=FALSE, sep='\t', row.names = FALSE)
					write.table(edge, file=paste(out[1], ".tsv", sep=""), quote=FALSE, sep='\t', row.names = FALSE)
					write.table(meta, file=paste(out[3], ".tsv", sep=""), quote=FALSE, sep='\t', row.names = FALSE)
				} else if(type %in% "xlsx") {
					stopifnot(length(out)==1)
					write.xlsx(node, file=paste(out, ".xlsx", sep=""), sheetName="Node Information", row.names=FALSE)
					write.xlsx(edge, file=paste(out, ".xlsx", sep=""), sheetName="Edge Information", row.names=FALSE)
					write.xlsx(meta, file=paste(out, ".xlsx", sep=""), sheetName="Metadata Information", row.names=FALSE)
				}
			} else {
				plot_data <- object@plot_data
				calc_type <- rep(object@calc_type, each=length(plot_data[,1]))
				calc_info <- rep(object@calc_info, each=length(plot_data[,1]))
				plot_data <- cbind(plot_data, calc_info, calc_type)
				if(type %in% "csv") {
					write.csv(plot_data, paste(out, ".csv", sep=""), row.names=FALSE)
				} else if(type %in% "tsv") {
					write.table(plot_data, file=paste(out, ".tsv", sep=""), quote=FALSE, sep='\t', row.names = FALSE)
				} else if(type %in% "xlsx") {
					write.xlsx(plot_data, file=paste(out, ".xlsx", sep=""), row.names=FALSE)
				}
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
#' @param type which type of output file. Can be "csv", "tsv" and "xlsx"
#' @return saves the data in the working directory as a csv and returns nothing
#' @export
setGeneric("save_analyser_object", function(object, which_data) standardGeneric("save_analyser_object"))
setMethod("save_analyser_object", "metime_analyser", function(object, which_data) {
			stopifnot(type %in% c("csv", "tsv", "xlsx"))
			list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
			list_of_col_data <- object@list_of_col_data[names(object@list_of_col_data) %in% which_data]
			list_of_row_data <- object@list_of_row_data[names(object@list_of_row_data) %in% which_data]
			for(i in 1:length(which_data)) {
				if(type %in% "csv") {
					write.csv(list_of_data[[i]], paste(which_data[i], ".csv", sep=""))
					write.csv(list_of_col_data[[i]], paste(which_data[i], "_col.csv", sep=""))
					write.csv(list_of_row_data[[i]], paste(which_data[i], "_row.csv", sep=""))	
				} else if(type %in% "tsv") {
					write.table(list_of_data[[i]], file=paste(which_data[i], ".tsv", sep=""), quote=FALSE, sep='\t', row.names = FALSE)
					write.table(list_of_col_data[[i]], file=paste(which_data[i], "_col.tsv", sep=""), quote=FALSE, sep='\t', row.names = FALSE)
					write.table(list_of_row_data[[i]], file=paste(which_data[i], "_row.tsv", sep=""), quote=FALSE, sep='\t', row.names = FALSE)
				} else if(type %in% "xlsx") {
					write.xlsx(list_of_data[[i]], file=paste(which_data[i], ".xlsx", sep=""), sheetName="data", row.names = FALSE)
					write.xlsx(list_of_col_data[[i]], file=paste(which_data[i], ".xlsx", sep=""), sheetName="coldata", row.names = FALSE)
					write.xlsx(list_of_row_data[[i]], file=paste(which_data[i], ".xlsx", sep=""), sheetName="rowdata", row.names = FALSE)
				}
					
			}
	}) 
