#' Get row_data from a S4 object of the class "metime_analyser" 
#' @description Get the row information of a specified dataset from a S4 object of the class "metime_analyser" 
#' @param object a S4 object of the class "metime_analyser" 
#' @param which_data a character of the dataset name of interest.
#' @return A dataframe containing the row data of the dataset of interest.
#' @export
setGeneric("get_rowdata", function(object, which_data) standardGeneric("get_rowdata"))
setMethod("get_rowdata", "metime_analyser", function(object, which_data) {
  if(!which_data %in% names(object@list_of_row_data)){
    warning(paste0("get_rowdata() could not find the data.frame ", which_data, "list_of_row_data of the object in the object"))
  }else{
    out <- object@list_of_row_data[[which_data]] %>% as.data.frame()
    return(out)
  }
	})


