
#' Modification (mod) function that merges a data frame to a dataset within the S4 object of class metime_analyser
#' @description a modification function that merges two sets of data including data, col_data and row_data.
#' @param object a S4 object of class metime_analyser.
#' @param which_data a vector of character defining which data should be merged. Has to contain two or more values.
#' @param type a character defining the data which should be appended. Can only be data, col_data, row_data. Default set to data
#' @param x a data.frame that is merged to the data.
#' @param id a character defining the column to be used for matching. Default set to NULL (first column of x will be used as id column)
#' @return a new S4 object of class metime_analyser with the  new merged dataset appended to it
#' @export 
#' 
setGeneric("add_data", function(object, which_data, type, x, id) standardGeneric("add_data")) 
setMethod("add_data", "metime_analyser", function(object, which_data, type="data", x, id=NULL)	{
  stopifnot(is.data.frame(x))
  stopifnot(type %in% c("data","col_data","row_data"))
  if(is.null(id)) id <- colnames(x)[1]
  
  if(type=="data"){
    object@list_of_data[[which_data]] <- object@list_of_data[[which_data]] %>% 
      dplyr::full_join(y = x, by=id) 
    for(i in setdiff(names(x), names(object@list_of_col_data[[which_data]]))){
      object@list_of_col_data[[which_data]][,i] <- NA
    }
  }else if(type=="col_data"){
    object@list_of_col_data[[which_data]] <- object@list_of_col_data[[which_data]] %>% 
      dplyr::full_join(y = x, by=id) 
  }else if(type=="row_data"){
    object@list_of_row_data[[which_data]] <- object@list_of_row_data[[which_data]] %>% 
      dplyr::full_join(y = x, by=id) 
  }
  # append information 
  out <- add_function_info(object=object, function_name="add_data", 
    params=list(which_data=which_data, x=colnames(x), id = id, type=type))
  
  return(out)
})
