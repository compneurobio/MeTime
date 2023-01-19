
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
#' @return data.frame with pairwise results
#' @export
setGeneric("calc_distance_pairwise", function(object, which_data, method) standardGeneric("calc_distance_pairwise"))
setMethod("calc_distance_pairwise", "metime_analyser", function(object, which_data, method="euclidean") {
  stopifnot(all(which_data %in% names(object@list_of_data)))
  
  flattenCorrMatrix <- function(cormat) {
    ut <- upper.tri(cormat)
    return(data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      dist  =(cormat)[ut]
    ))
  }
  
  my_data <-  lapply(which_data, function(x) object@list_of_data[[x]] %>% 
                       dplyr::mutate(id=rownames(.[]))) %>% 
    plyr::join_all(by="id", type="inner") %>% 
    `rownames<-`(.[,"id"]) %>% 
    dplyr::select(-id)
  
  if(dist %in% c("euclidean","maximum","minimum","manhattan","canberra","minkowski")){
    
    out <- my_data %>%
      stats::dist(method = method) %>%
      as.matrix() %>% 
      as.data.frame() %>% 
      flattenCorrMatrix() %>% 
      dplyr::mutate(type=method) %>%
      get_make_plotter_object(metadata=NULL, calc_type="pairwise_distance", 
                      calc_info = paste(which_data, "and" , method, "pairwise_distance", sep=" "),
                      plot_type="heatmap", style="ggplot", aesthetics=list(x="row", y="column", fill="dist"))
  }
  else {
    out=NA
  }
  return(out)
})


