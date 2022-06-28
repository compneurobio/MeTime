

#' Function to pack all the data into a single object of class "metab_analyser" 
#'
#' This function loads all the files from the parent directory. It assumes a 
#' certain naming pattern as follows: "datatype_[NULL|col|row]_data.rds" 
#' Any other naming pattern is not allowed. The function first writes 
#' all files into a list and each type of data is packed into its respective 
#' class i.e. col_data, row_data or phenotype data
#'
#' @param path Path to the parent directory
#' @return An object of class metab_analyser
#' @export

get_files_and_names <- function(path) {
	#path <- input$files$datapath
	data_list <- lapply(path, function(x) {
				x <- readRDS(x)
				return(x)
		})
	col_data_index <- grep("*_col_*", path)
	row_data_index <- grep("*_row_*", path)
	phenotype_index <- grep("*phenotype*", path)
	phenotype_row_index <- intersect(phenotype_index, row_data_index)
	phenotype_col_index <- intersect(phenotype_index, col_data_index)
	col_data_index <- col_data_index[!(col_data_index %in% phenotype_col_index)]
	row_data_index <- row_data_index[!(row_data_index %in% phenotype_row_index)]
	names_list <- lapply(path, function(x) {
					dummy <- unlist(lapply(strsplit(x, split="/"), function(b) return(b[length(b)])))
					dummy <- unlist(lapply(strsplit(dummy, split=".", fixed=TRUE), function(a) return(a[1])))
					return(dummy)
		})	
	list_of_data <- data_list[-c(col_data_index, row_data_index, phenotype_index)]
	list_of_col_data <- data_list[col_data_index]
	list_of_row_data <- data_list[row_data_index]
	names(list_of_data) <- names_list[-c(col_data_index, row_data_index, phenotype_index)]
	names(list_of_col_data) <- names_list[col_data_index]
	names(list_of_row_data) <- names_list[row_data_index]
	phenotype_index <- phenotype_index[!(phenotype_index %in% c(phenotype_row_index, phenotype_col_index))]
	phenotype <- as.data.frame(data_list[phenotype_index])
	metab_object <- new("metab_analyser", list_of_data=list_of_data, list_of_col_data=list_of_col_data, 
										list_of_row_data=list_of_row_data,
										phenotype=phenotype)
	return(metab_object)
}

#creating reference metab-analyser class that creates an object with full data

#' Constructor to generate an object of class metab_analyser. 
#' contains slots - list_of_data: For the list of all data matrices.
#'				  - list_of_col_data: list of all the col data files in the same order.
#' 				  - list_of_row_data: list of all the row data files in the same order.
#' 				  - phenotype: data.frame with phenotype data of the samples in it.  	
#' 	
#' @rdname metab_analyser
#' @export 
setClass("metab_analyser", slots=list(list_of_data="list", list_of_col_data="list", list_of_row_data="list", 
								 phenotype="data.frame")) 
