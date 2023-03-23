#' This function appends an object of class metime_analyser with a new dataset.
#' @description function to apply on metime_analyser object to append a new dataset into the existing object
#' @examples # append data frames into the metime_analyser object
#' appended_object <- add_dataset(object=metime_analyser_object, 
#'                    data=data, row_data=row_data, col_data=col_data, name="name of the new dataset")
#' @param object S4 object of class metime_analyser
#' @param data data.frame containing data 
#' @param col_data data.frame containing col_data: id column of col data has to match colnames of data
#' @param row_data data.frame containing row_data: id column of row data has to match rownames of data
#' @param name Name of the new dataset
#' @return An object of class metime_analyser
#' @export
setGeneric("add_dataset", function(object, data, col_data, row_data, name) standardGeneric("add_dataset"))
setMethod("add_dataset", "metime_analyser",function(object, data, col_data, row_data, name=NULL) {
  if(is.null(name)) name <- "set1"
  if(!all(rownames(data) %in% row_data$id) & !all(colnames(data) %in% col_data$id)) stop("id of col or row data do not match dataframe")
  if(!all(c("id","subject","time") %in% names(row_data))) stop("id, subject or timepoint column missing")

  if(!("covariates" %in% colnames(col_data))) {
  		col_data$covariates <- rep(NA, each=length(col_data$id))
  }
  
  object@list_of_data[[name]] <- data
  object@list_of_col_data[[name]] <- col_data
  object@list_of_row_data[[name]] <- row_data
  
  out <- object
  return(out)
})

