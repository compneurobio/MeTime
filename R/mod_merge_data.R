
#' Modification (mod) function to merge two sets of data or partial data merging for any analysis
#' @description a modification function that merges two sets of data or partial data including data, col_data and row_data.
#' If you want to use certain columns for the rest of the analysis properly for annotating plots then use mod_mutate 
#' function to rename the colnames and rownames (if applicable) to maintain consistency in the merged dataset 
#' @param object an S4 object of class metime_analyser.
#' @param which_data a vector of character defining which data should be merged. Has to contain two or more values.
#' If you want to merge partial data from other datasets then the first dataset will be consided to be complete unless and 
#' until you specify the first dataset as name in the cols_for_meta list. See below for example
#' @param filter_samples a character specifying if samples should be filtered. Default set to no filtering. 
#' Other options include common samples (samples used if found in all datasets), or set 
#' to name of dataset (samples filtered by samples of one dataset)
#' @param append logical. TRUE adds this new dataset to the object. FALSE creates a new object with 
#' just the data and the latest results section. Default is set to FALSE
#' @param cols_list A list of named character vectors to merge dataset. If the first dataset's name is used in the list
#' then those columns will be subsetted else the order of which_data will be followed irrespective of the names
#' example: which_data = c("dataset1", "dataset2"); cols_list = list(dataset1=c(...), dataset2=c(...)) or list(c(...))
#' If you want full data set cols_list=NULL
#' @param name a character vector to define the name of a new dataset.
#' @return a new S4 object of class metime_analyser with the  new merged dataset appended to it
#' @seealso [mod_merge_results]
#' @export 
setGeneric("mod_merge_data", function(object, which_data, filter_samples=NULL, cols_list=NULL, append=FALSE, name="merged_data") standardGeneric("mod_merge_data")) 
setMethod("mod_merge_data", "metime_analyser", function(object, which_data, filter_samples=NULL, cols_list=NULL, append=FALSE, name="merged_data")	{
  if(!length(which_data)>=2) {
    warning("At least two datasets are needed for merging. Exiting without making changes")
    return(object)
  }
  if(name %in% names(object@list_of_data)) {
    warning("Name has to be unique. Setting generic name for the time being")
    name <- "merged_data_1"
    if(name %in% names(object@list_of_data)) {
      index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
      index <- c(0:9)[grep(index, 0:9)+1]
      name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
    }
  }

  if(!is.null(names(cols_list)) & !is.null(cols_list)) {
    if(!which_data[1] %in% names(cols_list)) {
      message("dataset1's name is not there is cols_list checking order to proceed")
      if(!length(cols_list)==length(which_data)-1) {
        warning("Length doesn't match as well. Exiting without making any changes")
        return(object)
      }
      message("All okay moving on with merging")
    }
  } else if(!is.null(cols_list)) {
    if(!length(cols_list)==length(which_data)-1) {
      warning("Length of cols_list is not as expected. Exiting without making any changes")
      return(object)
    } else {
      names(cols_list) <- which_data
    }
  } else {
    cols_list <- list()
  }
  
  # data merge
  data <- lapply(which_data, function(x) {
    if(!x %in% names(cols_list)) {
      object@list_of_data[[x]] %>%
        dplyr::select(any_of(setdiff(names(object@list_of_data[[x]]), c("id","subject","time")))) %>% 
        dplyr::mutate(id = rownames(.[])) %>% return()
      } else {
        object@list_of_data[[x]] %>%
          dplyr::select(any_of(cols_list[[x]])) %>%
          dplyr::mutate(id=rownames(.[])) %>% return()
      }
    }) %>% 
    plyr::join_all(type="full", by="id") %>% 
    tibble::column_to_rownames(var="id")
  
  # row_data merge
  row_data <- lapply(which_data, function(x) {
      object@list_of_row_data[[x]] %>% return()
  }) %>% 
    plyr::join_all(by=c("id"),type="full") %>% 
    dplyr::filter(id %in% rownames(data))
  rownames(row_data) <- row_data$id
  
  # col_data merge
  col_data <- lapply(which_data, function(x) {
    if(!x %in% names(cols_list)) {
      object@list_of_col_data[[x]] %>% return()
     } else {
        object@list_of_col_data[[x]] %>% 
          dplyr::filter(id %in% cols_list[[x]]) %>% return()
     }
  }) %>%
    plyr::rbind.fill() %>% 
    dplyr::filter(id %in% colnames(data))
  rownames(col_data) <- col_data$id
  
  # filter samples
  if(!is.null(filter_samples)){
    if(filter_samples=="common"){
      use_rows <- lapply(which_data, function(x) object@list_of_row_data[[x]]$id) %>% 
        Reduce(f=intersect) %>% sort()
    } else if(filter_samples %in% names(object@list_of_row_data)) {
      use_rows = object@list_of_row_data[[filter_samples]]$id
    }
    data <- data[use_rows,]
    row_data <- row_data[which(row_data$id %in% use_rows),]
  }
  
  # append data to analyzer object
  if(append) {
		out <- add_dataset(object=object, data=data, col_data=col_data, row_data=row_data, name=name)
  } else {
    out <- get_make_analyser_object(data=data, col_data=col_data, row_data=row_data, name=name, 
      results=list(object@results[[length(object@results)]]))
  }
		
	# append information 
		out <- add_function_info(object=out, function_name="mod_merge_data", params=list(which_data=which_data,
      filter_samples=filter_samples, mutations=cols_list, name=name))
		return(out)
})
