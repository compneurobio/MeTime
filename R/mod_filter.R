
#' Modification function to filter columns in data, row_data or col_data
#' @description Modification (mod) function to filter columns in data, row_data or col_data. It is a wrapper function to dplyr::filter()
#' made for metime_analyser objects.
#' @param object A S4 object of class metime_analyser.
#' @param which_data a vector of character defining which data/result is to be filtered. 
#' @param type character input of length 1 to define the type of data to be manipulated. Accepted inputs are "row_data", "col_data", "data"
#' and "results". However renamed results will be returned to the user as a list of results and will not return the full analyser
#' object
#' @param ... arguments to pass directly into dplyr::filter() function.
#' @seealso [mod_mutate], [mod_rename]
#' @returns object with mutated data, col_data or row_data. However, if type is "results" then plot_data with modifications expected
#' will be returned
#' @export
setGeneric("mod_filter", function(object, which_data, ..., type="data") standardGeneric("mod_filter"))
setMethod("mod_filter", "metime_analyser", function(object, which_data, ..., type="data") {
    if(length(which_data)!=1) {
        warning("mod_filter(): length of which_data is not 1, exiting without making any changes")
        return(object)
    }
    if(!type %in% c("data", "row_data", "col_data", "results")) {
        warning("mod_filter(): type is unknown. Please check the accepted values and exiting without making any changes")
        return(object)
    }

    filter_exprs <- enquos(...)
    filter_exprs_str <- purrr::map_chr(filter_exprs, ~as.character(rlang::quo_text(.)))

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
    } else if(type %in% "results") {
        #if(is.null(results_index)) {
        #  warning("Results index is not specified, exiting without making any changes")
        #  return(object)
        #}
        if(grep("ggm|network", object@results[[which_data]]$information$calc_type) %>% length() >=1) {
          warning("mutations and filters are not possible for networks. Exiting without making any changes")
          return(object)
        }
        results <- object@results[[which_data]]
        results$plot_data <- lapply(seq_along(results$plot_data), function(x) {
            results$plot_data[[x]] <- results$plot_data[[x]] %>% dplyr::filter(!!!filter_exprs)
            return(results$plot_data[[x]])
        })
      return(results$plot_data)
    }
    
    object <- object %>% add_function_info(function_name="mod_filter", 
      params=list(which_data=which_data, type=type, mutations=paste(filter_exprs_str, collapse = ", ")))
    return(object)
  })