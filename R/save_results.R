#' Function to extract results into different types of files
#' @description extracts results from analyser object and saves it
#' @examples see examples here and maintain the order
#' If type is xlsx one out file is enough and the network data will be stored in different sheets
#' Network : save_results(object, results_index, out=c("edge", "node", "meta"), type="csv")
#' Others : save_results(object, results_index, out=c("outfile1", ...), type="tsv")
#' @param object An object of class metime_analyser
#' @param results_index character or numeric to define the results of interest
#' @param type character to define outfile type that is "csv", "xlsx" or "tsv"
#' @return saves the data into a csv and returns nothing
#' @export
setGeneric("save_results", function(object, results_index, type) standardGeneric("save_results"))
setMethod("save_results", "metime_analyser", function(object, results_index, type) {
			stopifnot(type %in% c("csv", "xlsx", "tsv"))
			results <- obejct@results[[results_index]]
			if(all(names(plot_data) %in% c("node", "edge", "metadata"))) {
				node <- results$plot_data$node 
				edge <-	results$plot_data$edge 
				meta <-	results$plot_data$metadata
				if(type %in% "csv") {
					write.csv(node, paste(results$information$calc_info, "_node.csv", sep=""), row.names=FALSE)
					write.csv(edge, paste(results$information$calc_info, "_edge.csv", sep=""), row.names=FALSE)
					write.csv(meta, paste(results$information$calc_info, "_metadta.csv", sep=""), row.names=FALSE)
				} else if(type %in% "tsv") {
					write.table(node, file=paste(results$information$calc_info, "_node.tsv", sep=""), quote=FALSE, sep='\t', row.names = FALSE)
					write.table(edge, file=paste(results$information$calc_info, "_node.tsv", sep=""), quote=FALSE, sep='\t', row.names = FALSE)
					write.table(meta, file=paste(results$information$calc_info, "_node.tsv", sep=""), quote=FALSE, sep='\t', row.names = FALSE)
				} else if(type %in% "xlsx") {
					xlsx::write.xlsx(node, file=paste(results$information$calc_info, ".xlsx", sep=""), sheetName="Node Information", row.names=FALSE)
					xlsx::write.xlsx(edge, file=paste(results$information$calc_info, ".xlsx", sep=""), sheetName="Edge Information", row.names=FALSE)
					xlsx::write.xlsx(meta, file=paste(results$information$calc_info, ".xlsx", sep=""), sheetName="Metadata Information", row.names=FALSE)
				}
			} else {
				lapply(seq_along(results$plot_data), function(x) {
						if(type %in% "csv") {
							write.csv(results$plot_data[[x]], 
								paste(results$information$calc_info[x], ".csv", sep=""), row.names=FALSE)
						} else if(type %in% "tsv") {
							write.table(results$plot_data[[x]], 
								file=paste(results$information$calc_info[x], ".tsv", sep=""), 
								quote=FALSE, sep='\t', row.names = FALSE)
						} else if(type %in% "xlsx") {
							xlsx::write.xlsx(results$plot_data[[x]], 
								file=paste(".xlsx", sep=""), row.names=FALSE)
						}
					})
			}
	}) 


