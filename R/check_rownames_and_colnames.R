

#' Function to check the format of rownames and colnames and if they are same or not
#' @description sanity check to check for rownames of the data
#' @param object S4 object of class of metime_analyser
#' @return NULL if it passes all the sanity checks
#' @export
setGeneric("check_rownames_and_colnames", function(object) standardGeneric("check_rownames_and_colnames"))
setMethod("check_rownames_and_colnames", "metime_analyser", function(object) {
		list_of_data <- object@list_of_data
		list_of_row_data <- object@list_of_row_data
		list_of_col_data <- object@list_of_col_data
		#data_checks for duplicates
		for(i in 1:length(list_of_data)) {
			if(length(colnames(list_of_data[[i]])) != length(unique(colnames(list_of_data[[i]])))) {
				stop(paste("In Dataset", names(list_of_data)[i], "some metabolites are duplicated", sep=" "))
			} else {
				print(paste("In Dataset", names(list_of_data)[i], "NO metabolites are duplicated and dataset is good", sep=" "))
			} 
			if(length(rownames(list_of_data[[i]])) != length(unique(rownames(list_of_data[[i]])))) {
				stop(paste("In Dataset", names(list_of_data)[i], "some samples are duplicated", sep=" "))	
			} else {
				print(paste("In Dataset", names(list_of_data)[i], "NO samples are duplicated and dataset is good", sep=" "))
			}
		}
		#data checks for colnames and rownames in col_data and row_data respectively
		for(i in 1:length(list_of_data)) {
			if(!(names(list_of_data)[i] %in% object@annotations$phenotype | names(list_of_data)[i] %in% object@annotations$medication)) {
				if(all(list_of_col_data[[i]]$id %in% colnames(list_of_data[[i]]))) {
					print(paste("All metablites are matched in col_data and columns of raw data in dataset", names(list_of_data)[i], sep=" "))
				} else {
					stop(paste("In Dataset", names(list_of_data)[i], "all metabolites are not mapped in the col data"))
				}
				if(all(list_of_row_data[[i]]$id %in% rownames(list_of_data[[i]]))) {
					print(paste("All samples are matched in row_data and rows of raw data in dataset", names(list_of_data)[i], sep=" "))
				} else {
					stop(paste("In Dataset", names(list_of_data)[i], "all samples are not mapped in the row data"))
				}
			}
		}	
	})


