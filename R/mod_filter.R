#' Modification function to filter columns in data, row_data or col_data
#' @description Modification (mod) function to filter columns in data, row_data or col_data
#' @param object A S4 object of class metime_analyser.
#' @param which_data a vector of character defining which data should be merged. 
#' @param type a character defining which data to filter. Can be "data","row_data", or "col_data". Default set to "data". Note: if a dataset is filtered the row_data will be changed accordingly
#' @param ... arguments to pass directly into dplyr::filter() function.
#' @returns object with mutated data, col_data or row_data
#' @export
setGeneric("mod_filter", function(object, which_data, type="data", ...) standardGeneric("mod_filter"))
setMethod("mod_filter", "metime_analyser", function(object, which_data, type="data", ...) {
    stopifnot(length(which_data)==1)
    stopifnot(which_data %in% names(object@list_of_data))
    stopifnot(type %in% c("data", "row_data", "col_data"))

    filter_exprs <- enquos(...)
    filter_exprs_str <- purrr::map_chr(filter_exprs, ~as.character(quo_text(.)))

    if(type %in% "data") {
      #filter dataset
        data <- object %>% get_data(which_data=which_data)
        data <- data %>% dplyr::filter(!!!filter_exprs)
        object@list_of_data[[which_data]] <- data
        #adjust row data accordingly
        object@list_of_row_data[[which_data]] <- object %>% 
        get_rowdata(which_data=which_data) %>%  
        dplyr::filter(id %in% rownames(object@list_of_data[[which_data]]))  
    } else if(type %in% "col_data") {
        data <- object %>% get_coldata(which_data=which_data)
        data <- data %>% dplyr::filter(!!!filter_exprs)
        object@list_of_col_data[[which_data]] <- data
        #adjust data accordingly
        object@list_of_data[[which_data]] <- object %>% 
        get_data(which_data=which_data) %>% 
        dplyr::select(any_of(object@list_of_col_data[[which_data]]$id))
    } else if(type %in% "row_data") {
        data <- object %>% get_rowdata(which_data=which_data)
        data <- data %>% dplyr::filter(!!!filter_exprs)
        object@list_of_row_data[[which_data]] <- data
        #adjust data accordingly
        object@list_of_data[[which_data]] <- object %>% 
        get_data(which_data=which_data) %>%
        dplyr::mutate(id = rownames(.[])) %>% 
        dplyr::filter(id %in% object@list_of_row_data[[which_data]]$id) %>% 
        dplyr::select(-id)
    }
    
    object <- object %>% add_function_info(function_name="mod_filter", 
      params=list(which_data=which_data, type=type, mutations=paste(filter_exprs_str, collapse = ", ")))
    return(object)
  })