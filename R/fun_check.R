

#' Function to check the ids in the data and data format
#' @description sanity check to check for ids and order of the data
#' @param object S4 object of class of metime_analyser
#' @return NULL if it passes all the sanity checks
#' @export
setGeneric("check_ids_and_classes",  function(object) standardGeneric("check_ids_and_classes")) 
setMethod("check_ids_and_classes", "metime_analyser", function(object) {
			list_of_interest <- list(data=object@list_of_data, col_data=object@list_of_col_data, row_data=object@list_of_row_data)
			names <- names(object@list_of_data)
			type <- unlist(lapply(strsplit(names(object), split="list_of_"), function(x) return(x[2])))
			for(i in 1:length(list_of_interest)) {
				for(j in 1:length(list_of_interest[[i]])) {
					if(class(!list_of_interest[[i]][[j]]) %in% "data.frame") {
						stop(paste("In dataset", names[j], "datatype", type[i], 
							"class of data should be data.frame only"))
					} 
					if(!all(rownames(list_of_interest[[i]][[j]] == list_of_interest[[i]][[j]]$id))) {
						stop(paste("In dataset", names[j], "datatype", type[i], 
							"rownames and ids are not mapped correctly"))
					}
				}
			}
			return(NULL)
	})


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

#' Function to check for col_normality data whether it is added or not. 
#' @description function to check whether col_normality data is added to
#' the object or not
#' @param object S4 object of class of metime_analyser
#' @param which_data dataset/s to check
#' @return NULL if it passes all the sanity checks
#' @export
setGeneric("check_col_normality", function(object, which_data) standardGeneric("check_col_normality"))
setMethod("check_col_normality", "metime_analyser", function(object, which_data) {
			check_data <- object@list_of_col_data[names(object@list_of_col_data) %in% which_data]
			cols <- c("pval", "statistic", "normal")
			for(i in 1:length(check_data)) {
					if(is.null(grep(paste(cols, collapse="|"), colnames(check_data[[i]])))) {
						print(paste("normality data is not added to", which_data[i], "dataset", sep=" "))
					} else {
						if(is.null(grep("shapiro", colnames(check_data[[i]])))) {
							print(paste("shapiro normality test is not added to", which_data[i], "dataset", sep=" "))
						} else if(is.null(grep("kruskal", colnames(check_data[[i]])))) {
							print(paste("kruskal normality test is not added to", which_data[i], "dataset", sep=" "))
						} else {
							print(paste("All normality tests are added to", which_data[i], "dataset", sep=" "))
						}
					}
			}
			return(NULL)
	})