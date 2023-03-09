#' Function to extract analyser object data into a csv
#' @description extracts information from analyser object and saves it as a csv
#' @examples see examples here
#' save_analyser_object(object, which_data="dataset")
#' @param object An object of class metime_plotter
#' @param which_data Character to specify the dataset
#' @param type which type of output file. Can be "csv", "tsv" and "xlsx"
#' @return saves the data in the working directory as a csv and returns nothing
#' @export
setGeneric("save_analyser_object", function(object, which_data, type) standardGeneric("save_analyser_object"))
setMethod("save_analyser_object", "metime_analyser", function(object, which_data, type) {
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
					xlsx::write.xlsx(list_of_data[[i]], file=paste(which_data[i], ".xlsx", sep=""), sheetName="data", row.names = FALSE)
					xlsx::write.xlsx(list_of_col_data[[i]], file=paste(which_data[i], ".xlsx", sep=""), sheetName="coldata", row.names = FALSE)
					xlsx::write.xlsx(list_of_row_data[[i]], file=paste(which_data[i], ".xlsx", sep=""), sheetName="rowdata", row.names = FALSE)
				}
					
			}
	}) 
