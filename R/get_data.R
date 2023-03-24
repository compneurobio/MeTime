#' Get data from a S4 object of the class "metime_analyser" 
#' @description Get the data of a specified dataset from a S4 object of the class "metime_analyser" 
#' @param object a S4 object of the class "metime_analyser" 
#' @param which_data a character of the dataset name of interest.
#' @return A dataframe containing the data of the dataset of interest.
#' @export
setGeneric("get_data", function(object, which_data) standardGeneric("get_data"))
setMethod("get_data", "metime_analyser", function(object, which_data) {
  if(!which_data %in% names(object@list_of_data)){
    warning(paste0("get_data() could not find the data.frame ", which_data, "list_of_data of the object"))
  }else{
    out <- object@list_of_data[[which_data]] %>% as.data.frame()
    return(out)
  }
	})