
#' Function to calculate dissimilarity using distance measures 
#'
#' @description calculate pairwise distances between samples 
#' This function creates a dataframe for plotting from a dataset.
#' @examples # Example to calculate pairwise distances
#' dist <- calc_pairwise_distance(object=metime_analyser_object, which_data="name of the dataset", 
#'           method="euclidean")
#' @param object S4 Object of class metime_analyser
#' @param which_data specify datasets to calculate on. One or more possible
#' @param method default setting: method="euclidean", Alternative "maximum","minimum",
#' "manhattan","canberra","minkowski" are also possible
#' @param name name of the results should be of length=1
#' @param stratifications List to stratify data into a subset. Usage list(name=value)
#' @return data.frame with pairwise results
#' @export
setGeneric("calc_distance_samples", function(object, which_data, method, name="calc_distance_pairwise_1", stratifications) standardGeneric("calc_distance_samples"))
setMethod("calc_distance_samples", "metime_analyser", function(object, which_data, method="euclidean", name="calc_distance_pairwise_1", stratifications) {
  stopifnot(all(which_data %in% names(object@list_of_data)))
  if(grep(name, names(object@results)) %>% length() >=1) {
    warning("name of the results was previously used, using a different name")
    index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
    index <- c(0:9)[grep(index, 0:9)+1]
    name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
  }
  
  flattenCorrMatrix <- function(cormat) {
    ut <- upper.tri(cormat)
    return(data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      dist  =(cormat)[ut]
    ))
  }
  if(!is.null(stratifications) && length(stratifications) >= 1) {
    data_list <- get_stratified_data(object=object, which_data=which_data, stratifications=stratifications)
    my_data <- data_list[["data"]]
  } else {
    my_data <-  lapply(which_data, function(x) object@list_of_data[[x]] %>% 
                         dplyr::mutate(id=rownames(.[]))) %>% 
      plyr::join_all(by="id", type="inner") %>% 
      `rownames<-`(.[,"id"]) %>% 
      dplyr::select(-id)
  }
  if (nrow(my_data) < 2 || ncol(my_data) < 2) {
    warning("calc_distance_samples(): need at least 2 samples and 2 features after filtering")
    return(object)
  }
  
  if(method %in% c("euclidean","maximum","minimum","manhattan","canberra","minkowski")){
    
    out <- my_data %>%
      stats::dist(method = method) %>%
      as.matrix() %>%   
      as.data.frame() %>% 
      flattenCorrMatrix() %>% 
      dplyr::mutate(type=method)
    out <- get_make_results(data=list(pairwise_distance=out), object=object, metadata=NULL, calc_type="pairwise_distance", 
                      calc_info = paste(which_data, "and" , method, "pairwise_distance_samples", sep=" "),
                      name=name) %>%
        add_function_info(function_name="calc_distance_samples",
                params=list(which_data=which_data, method=method, stratifications=stratifications))
  }
  else {
    out=NA
    warning("calc_distance_samples(): unsupported method")
    return(object)
  }
  return(out)
})


