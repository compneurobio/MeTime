
#' Function to calculate dissimilarity using distance measures 
#'
#' @description calculate pairwise distances
#' This function creates a dataframe for plotting from a dataset.
#' @examples # Example to calculate pairwise distances
#' dist <- calc_pairwise_distance(object=metime_analyser_object, which_data="name of the dataset", 
#'           method="euclidean")
#' @param object S4 Object of class metime_analyser
#' @param which_data specify datasets to calculate on. One or more possible
#' @param method default setting: method="euclidean", Alternative "maximum","minimum",
#' "manhattan","canberra","minkowski" are also possible
#' @param name name of the results should be of length=1
#' @param cols_for_meta list equal to length of which_data defining the columns for metadata
#' @param stratifications List to stratify data into a subset. Usage list(name=value)
#' @return data.frame with pairwise results
#' @export
setGeneric("calc_distance_pairwise", function(object, which_data, method, name, cols_for_meta, stratifications) standardGeneric("calc_distance_pairwise"))
setMethod("calc_distance_pairwise", "metime_analyser", function(object, which_data, method="euclidean", name, cols_for_meta, stratifications) {
  stopifnot(all(which_data %in% names(object@list_of_data)))
  
  flattenCorrMatrix <- function(cormat) {
    ut <- upper.tri(cormat)
    return(data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      dist  =(cormat)[ut]
    ))
  }
  if(is.null(cols_for_meta)) {
      metadata <- NULL
  } else {
       metadata <- get_metadata_for_columns(object = object, 
                                         which_data = which_data, 
                                         columns = cols_for_meta, 
                                         names = c("name", "group"), 
                                         index_of_names = "id")
  }
  
  my_data <-  lapply(which_data, function(x) object@list_of_data[[x]] %>% 
                       dplyr::mutate(id=rownames(.[]))) %>% 
    plyr::join_all(by="id", type="inner") %>% 
    `rownames<-`(.[,"id"]) %>% 
    dplyr::select(-id)
  if(length(stratifications)>=1) {
        dummy_data <- object@list_of_data[[which_data[1]]]
        row_data <- object@list_of_row_data[[which_data[1]]]
        stratifications <- lapply(names(stratifications), function(x) {
              row_data <- row_data[row_data[,x] %in% stratifications[[x]], ]
              return(stratifications[[x]]) 
          })
        my_data <- my_data[rownames(my_data) %in% rownames(row_data), ]
  }
  
  if(dist %in% c("euclidean","maximum","minimum","manhattan","canberra","minkowski")){
    
    out <- my_data %>%
      stats::dist(method = method) %>%
      as.matrix() %>% 
      as.data.frame() %>% 
      flattenCorrMatrix() %>% 
      dplyr::mutate(type=method)
    out <- get_make_results(data=list(out), object=object, metadata=metadata, calc_type="pairwise_distance", 
                      calc_info = paste(which_data, "and" , method, "pairwise_distance", sep=" "),
                      name=name)
  }
  else {
    out=NA
  }
  return(out)
})


