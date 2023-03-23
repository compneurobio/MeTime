#' Add a data.frame to analyzer object. 
#' @description Add a data.frame x to existing data within an analyzer object. Data can be added to data, col_data or row_data. A full join of the data is used.
#' @param object a S4 object of the class "metime_analyzer".
#' @param which_data a character to define which dataset is to be used.
#' @param type a character defining the data which should be appended. Can only be data, col_data, row_data. Default set to data.
#' @param x a data.frame that is merged to the data.
#' @param id a character vector defining the column to be used for matching. Default set to NULL (first column of x will be used as id column)
#' @return a S4 object of class "metime_analyser" with appended data.frame x
#' @export 
setGeneric("add_data", function(object, which_data, type, x, id) standardGeneric("add_data")) 
setMethod("add_data", "metime_analyser", function(object, which_data, type="data", x, id=NULL)	{
  out <- object 
  if(!is.data.frame(x)) warning("add_data() x is not a data.frame.")
  else if(!type %in% c("data","col_data","row_data")) warning("add_data() type can only be 'data', 'col_data' or 'row_data'")
  else{
    if(is.null(id)) id <- colnames(data)[1]
    # if id is null then the first column of the data.frame is used as id
    if(type=="data"){
    out@list_of_data[[which_data]] <- out@list_of_data[[which_data]] %>% 
      dplyr::full_join(y = x, by=id) 
    for(i in setdiff(names(x), names(out@list_of_col_data[[which_data]]))){
      out@list_of_col_data[[which_data]][,i] <- NA
    }
  }else if(type=="col_data"){
    out@list_of_col_data[[which_data]] <- out@list_of_col_data[[which_data]] %>% 
      dplyr::full_join(y = x, by=id) 
  }else if(type=="row_data"){
    out@list_of_row_data[[which_data]] <- out@list_of_row_data[[which_data]] %>% 
      dplyr::full_join(y = x, by=id) 
  }
  # append information 
  out <- add_function_info(object=out, function_name="add_data", 
    params=list(which_data=which_data, x=colnames(x), id = id, type=type))
  }
  return(out)
})
