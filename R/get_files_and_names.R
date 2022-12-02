#' Function to pack all the data into a single object of class "metime_analyser" 
#'
#' @description This function loads all the files from the parent directory. It assumes a 
#' certain naming pattern as follows: "datatype_None|col|row_data.rds" 
#' Any other naming pattern is not allowed. The function first writes 
#' all files into a list and each type of data is packed into its respective 
#' class i.e. col_data, row_data or data
#'
#' @examples

#' # Input in the parent directory from which the data files are to be extracted along with annotations_index to specify phenotype and medication data

#' get_files_and_names(path=/path/to/parent/directory, annotations_index=list(phenotype="Name of phenotype file", medication="name of phenotype file"))

#' @param path Path to the parent directory
#' @param annotations_index a list to be filled as follows = list(phenotype="Name or index of the files", medication="Name or index of the files")
#' @return An object of class metime_analyser
#' @export

get_files_and_names <- function(path, annotations_index) {
	#order the names and columns
	#Add timepoints and subjects to row data
	#check the order because order subjects first and the timepoints next
	#Add time column as "time" 
	#path <- input$files$datapath
	path <- list.files(path, pattern="[.rds|.RDS]", full.names=TRUE)	
	data_list <- lapply(path, function(x) {
				x <- readRDS(x)
				return(x)
		})
	col_data_index <- grep("*_col_*", path)
	row_data_index <- grep("*_row_*", path)
	names_list <- lapply(path, function(x) {
					dummy <- unlist(lapply(strsplit(x, split="/"), function(b) return(b[length(b)])))
					dummy <- unlist(lapply(strsplit(dummy, split=".", fixed=TRUE), function(a) return(a[1])))
					return(dummy)
		})	
	list_of_data <- data_list[-c(col_data_index, row_data_index)]
	list_of_col_data <- data_list[col_data_index]
	list_of_row_data <- data_list[row_data_index]
	names(list_of_data) <- names_list[-c(col_data_index, row_data_index)]
	names(list_of_col_data) <- names_list[col_data_index]
	names(list_of_col_data) <- gsub("_col_data", "_data", names(list_of_col_data))
	names(list_of_row_data) <- names_list[row_data_index]
	names(list_of_row_data) <- gsub("_row_data", "_data", names(list_of_row_data))
	metab_object <- new("metime_analyser", list_of_data=list_of_data, list_of_col_data=list_of_col_data, 
										list_of_row_data=list_of_row_data,
										annotations=annotations_index)
	#sanity checks for the object created
	check_rownames_and_colnames(metab_object)
	#Update subject and time columns in the row data
	for(i in 1:length(metab_object@list_of_row_data)) {
			if("rid" %in% colnames(metab_object@list_of_row_data[[i]])) {
				metab_object@list_of_row_data[[i]]$rid <- NULL
			}
			if("timepoint" %in% colnames(metab_object@list_of_row_data[[i]])) {
				metab_object@list_of_row_data[[i]]$timepoint <- NULL
			}
			if("RID" %in% colnames(metab_object@list_of_row_data[[i]])) {
				metab_object@list_of_row_data[[i]]$RID <- NULL
			}
	}
	for(i in 1:length(metab_object@list_of_row_data)) {
			if("subject" %in% colnames(metab_object@list_of_row_data[[i]]) && "time" %in% colnames(metab_object@list_of_row_data[[i]])) {
					metab_object@list_of_row_data[[i]] <- metab_object@list_of_row_data[[i]] %>% arrange(subject, time)
			} else {
					subject <- unlist(lapply(strsplit(rownames(metab_object@list_of_data[[i]]), split="_"), function(x) return(x[1])))
					time <- unlist(lapply(strsplit(rownames(metab_object@list_of_data[[i]]), split="_"), function(x) return(x[2])))
					metab_object@list_of_row_data[[i]] <- as.data.frame(cbind(metab_object@list_of_row_data[[i]], subject, time))
			}
			metab_object@list_of_data[[i]] <- metab_object@list_of_data[[i]][order(rownames(metab_object@list_of_row_data[[i]])), ]
	}
	for(i in 1:length(metab_object@list_of_col_data)) {
			metab_object@list_of_data[[i]] <- metab_object@list_of_data[[i]][ ,order(colnames(metab_object@list_of_data[[i]]))]
			metab_object@list_of_col_data[[i]] <- metab_object@list_of_col_data[[i]][order(metab_object@list_of_col_data[[i]]$id), ]
			if("covariates" %in% colnames(metab_object@list_of_col_data[[i]])) {
					next
			} else {
					metab_object@list_of_col_data[[i]]$covariates <- rep(NA, each=length(metab_object@list_of_col_data[[i]]$id))
			}
	}
	out <- metab_object
	return(out)
}
