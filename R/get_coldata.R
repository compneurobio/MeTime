#' Get col_data from a S4 object of the class "metime_analyser" 
#' @description Get the column information of a specified dataset from a S4 object of the class "metime_analyser" 
#' @param object a S4 object of the class "metime_analyser" 
#' @param which_data a character of the dataset name of interest.
#' @return A dataframe containing the col data of the dataset of interest.
#' @export
setGeneric("get_coldata", function(object, which_data) standardGeneric("get_coldata"))
setMethod("get_coldata", "metime_analyser", function(object, which_data) {
  if(!which_data %in% names(object@list_of_col_data)){
    warning(paste0("get_coldata() could not find the data.frame ", which_data, "list_of_col_data of the object"))
  }else{
    out <- object@list_of_col_data[[which_data]] %>% as.data.frame()
    return(out)
  }
	})