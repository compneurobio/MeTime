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
					xlsx::write.xlsx(node, file=paste(out, ".xlsx", sep=""), sheetName="Node Information", row.names=FALSE)
					xlsx::write.xlsx(edge, file=paste(out, ".xlsx", sep=""), sheetName="Edge Information", row.names=FALSE)
					xlsx::write.xlsx(meta, file=paste(out, ".xlsx", sep=""), sheetName="Metadata Information", row.names=FALSE)
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
					xlsx::write.xlsx(plot_data, file=paste(out, ".xlsx", sep=""), row.names=FALSE)
				}
			}
	}) 


