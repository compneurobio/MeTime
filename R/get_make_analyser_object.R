#' Function to pack all the data into a single object of class "metime_analyser" 
#'
#' @description This function creates an object of class metime_analyser from a dataset.
#' @examples
#' # new_metime_analyser_object <- get_make_metab_object(data=data_frame, col_data=col_data_frame, row_data=row_data, name="name of the new dataset", 
#'                                annotations_index=list(phenotype="name of phenotype", medication="name of medication"))
#' @param data data.frame containing data 
#' @param col_data data.frame containing col_data: id column of col data has to match colnames of data
#' @param row_data data.frame containing row_data: id column of row data has to match rownames of data
#' @param annotations_index a list to be filled as follows = list(phenotype="Name or index of the file/list", medication="Name or index of the files/list")
#' @param name character. Name you want to assign to the new dataset that is being added on
#' @return An object of class metime_analyser
#' @export

get_make_analyser_object <- function(data, col_data, row_data, annotations_index=list(), name=NULL) {
  if(is.null(name)) name <- "set1"
  if(!all(rownames(data) %in% row_data$id) & !all(colnames(data) %in% col_data$id)) stop("id of col or row data do not match dataframe")
  if(!all(c("id","subject","time") %in% names(row_data))) stop("id, subject or time column missing")
  
  list_of_data <- list()
  list_of_data[[name]] <- data
  
  list_of_col_data <- list()
  if(!("covariates" %in% colnames(col_data))) {
  		col_data$covariates <- rep(NA, each=length(col_data$id))
  }
  list_of_col_data[[name]] <- col_data 
  
  list_of_row_data <- list()
  list_of_row_data[[name]] <- row_data
  
  out <- new("metime_analyser", 
                      list_of_data=list_of_data, 
                      list_of_col_data=list_of_col_data, 
                      list_of_row_data=list_of_row_data,
                      annotations=annotations_index)
  
  return(out)
}

