
#' Modification (mod) function to merge two sets of data for any analysis
#' @description a modification function that merges two sets of data including data, col_data and row_data.
#' @param object a S4 object of class metime_analyser.
#' @param which_data a vector of character defining which data should be merged. Has to contain two or more values.
#' @param filter_samples a character specifying if samples should be filtered. Default set to no filtering. Other options include common samples (samples used if found in all datasets), or set to name of dataset (samples filtered by samples of one dataset)
#' @param name a character vector to define the name of a new dataset.
#' @importClassesFrom metime_analyser
#' @return a new S4 object of class metime_analyser with the  new merged dataset appended to it
#' @export 
setGeneric("mod_merge_data", function(object, which_data, filter_samples=NULL, name="merged_data") standardGeneric("mod_merge_data")) 
setMethod("mod_merge_data", "metime_analyser", function(object, which_data, filter_samples=NULL, name="merged_data")	{
  stopifnot(length(which_data)>=2) # two or more data index are needed for merging
  stopifnot(!name %in% names(object@list_of_data)) # name has to be unique
  
  # data merge
  data <- lapply(which_data, function(x)
    object@list_of_data[[x]] %>%
      dplyr::select(any_of(setdiff(names(object@list_of_data[[x]]), c("id","subject","time")))) %>% 
      dplyr::mutate(id = rownames(.[]))
    ) %>% 
    plyr::join_all(type="full", by="id") %>% 
    tibble::column_to_rownames(var="id")
  
  # row_data merge
  row_data <- lapply(which_data, function(x)
    object@list_of_row_data[[x]]
  ) %>% 
    plyr::join_all(by=c("id"),type="full") %>% 
    dplyr::filter(id %in% rownames(data))
  
  # col_data merge
  col_data <- lapply(which_data, function(x)
    object@list_of_col_data[[x]]
  ) %>%
    plyr::rbind.fill() %>% 
    dplyr::filter(id %in% colnames(data))
  
  # filter samples
  if(!is.null(filter_samples)){
    if(filter_samples=="common"){
      use_rows <- lapply(which_data, function(x) object@list_of_row_data[[x]]$id) %>% 
        Reduce(f=intersect) %>% sort()
    }else if(filter_samples %in% names(object@list_of_row_data)){
      use_rows = object@list_of_row_data[[filter_samples]]$id
    }
    data <- data[use_rows,]
    row_data <- row_data[which(row_data$id %in% use_rows),]
  }
  
  # append data to analyzer object
		out <- get_append_analyser_object(object=object, data=data, col_data=col_data, row_data=row_data, name=name)
		
	# append information 
		out <- add_function_info(object=out, function_name="mod_merge_data", params=list(which_data=which_data,filter_samples=filter_samples, name=name))
		return(out)
})
