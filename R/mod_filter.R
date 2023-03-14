
#' Modification function to filter columns in data, row_data or col_data
#' @description Modification (mod) function to filter columns in data, row_data or col_data
#' @param object A S4 object of class metime_analyser.
#' @param which_data a vector of character defining which data/result is to be filtered. 
#' @param type either "row" for row_data or "col" for col_data or "data" for data and "results" for results. 
#' Set to "data" as default
#' @param ... arguments to pass directly into dplyr::filter() function.
#' @returns object with mutated data, col_data or row_data
#' @export
setGeneric("mod_filter", function(object, which_data, ..., type="data") standardGeneric("mod_filter"))
setMethod("mod_filter", "metime_analyser", function(object, which_data, ..., type="data") {
    stopifnot(length(which_data)==1)
    stopifnot(which_data %in% names(object@list_of_data))
    stopifnot(type %in% c("data", "row_data", "col_data", "results"))


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