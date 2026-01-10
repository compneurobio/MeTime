#' Pack all the data into a single object of class "metime_analyser" 
#' @description Create an object of class "metime_analyser" from a dataset including data, row_data (id column corresponds to rownames(data)) and col_data (id colummn corresponds to colnames(data))
#' @param data a data.frame containing data.
#' @param col_data a data.frame containing col_data: id column of col data has to match colnames of data.
#' @param row_data a data.frame containing row_data: id column of row data has to match rownames of data.
#' @param annotations_index a named list to be filled as list(phenotype="Name or index of the file/list", medication="Name or index of the files/list").
#' @param name a character to be assign to the new dataset. Default is set to "set_1".
#' @param results a list of existing results. Default set to NULL.
#' @return An object of class metime_analyser
#' @export
get_make_analyser_object <- function(data, col_data, row_data, annotations_index=list(), name="set_1", results=list()) {
  if(!all(rownames(data) %in% row_data$id) | !all(colnames(data) %in% col_data$id)) stop("id of col or row data do not match dataframe")
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
  annotations_index <- list(annotations_index, modifications=list())
  
  out <- new("metime_analyser", 
                      list_of_data=list_of_data, 
                      list_of_col_data=list_of_col_data, 
                      list_of_row_data=list_of_row_data,
                      annotations=annotations_index,
                      results=results)
  return(out)
}
