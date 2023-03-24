#' Pack a dataset (data, row_data, col_data) into a single object of class "metime_analyser".
#' @description This function loads all the files from the parent directory. It assumes a 
#' certain naming pattern as follows: "datatype_None|col|row_data.rds" 
#' Any other naming pattern is not allowed. The function first writes 
#' all files into a list and each type of data is packed into its respective 
#' class i.e. col_data, row_data or data. 
#' @param path a character defining the path to the parent directory
#' @param annotations_index a named list to be filled as 
#' list(phenotype="Name or index of the files", 
#'	medication="Name or index of the files")
#' @examples
#' get_files_and_names(path="/path/to/parent/directory", 
#' 	annotations_index=list(phenotype="Name of phenotype file", 
#'	medication="name of phenotype file"))
#' @return An object of class "metime_analyser".
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
	annotations_index <- list(annotations_index, modifications=list())
	metab_object <- new("metime_analyser", list_of_data=list_of_data, list_of_col_data=list_of_col_data, 
										list_of_row_data=list_of_row_data,
										annotations=annotations_index,
										results=list())
	#Update subject and time columns in the row data
	out <- metab_object
	out@list_of_row_data <- lapply(seq_along(out@list_of_row_data), function(x) {
			a <- out@list_of_row_data[[x]]
			if("rid" %in% colnames(a)) {
				a$rid <- NULL
			} else if("timepoint" %in% colnames(a)) {
				a$timepoint <- NULL
			} else if("RID" %in% colnames(a)) {
				a$RID <- NULL
			}
			if("subject" %in% colnames(a) && "time" %in% colnames(a)) {
				a <- a %>% dplyr::arrange(subject, time)
			} else {
				a$subject <- rownames(a) %>% gsub(pattern="_[a-z|A-Z][0-9]+", replacement="")
				a$time <- rownames(a) %>% gsub(pattern="[a-z|A-Z][0-9]+_", replacement="")
			}
			a <- a[order(rownames(a)), ]
			return(a)
	})
	names(out@list_of_row_data) <- names(out@list_of_data)
	out@list_of_col_data <- lapply(seq_along(out@list_of_col_data),function(x) {
			a <- out@list_of_col_data[[x]]
			if(!"covariates" %in% colnames(a)) {
				a$covariates <- rep(NA, each=length(a$id))
			}
			return(a)
	})
	names(out@list_of_col_data) <- names(out@list_of_data)
	return(out)
}



