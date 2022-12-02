
#' Function to calculate correlation 
#'
#' @description calculate pairwise correlations
#' This function creates a dataframe for plotting from a dataset.
#' @examples # Example to calculate correlations
#' dist <- calc_correlation(object=metime_analyser_object, which_data="name of the dataset", 
#'           method="pearson")
#' @param object S4 Object of class metime_analyser
#' @param which_data specify datasets to calculate on. One or more possible
#' @param method default setting: method="pearson", Alternative "spearman" also possible
#' @return data.frame with pairwise results
#' @export
setGeneric("calc_correlation_pairwise", function(object, which_data, method) standardGeneric("calc_correlation_pairwise"))
setMethod("calc_correlation_pairwise", "metime_analyser", function(object, which_data, method="pearson"){
  stopifnot(all(which_data %in% names(object@list_of_data)))
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    return(data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    ))
  }
    
  my_data <-  lapply(which_data, function(x) object@list_of_data[[x]] %>% 
                         dplyr::mutate(id=rownames(.[]))) %>% 
      plyr::join_all(by="id", type="inner") %>% 
      `rownames<-`(.[,"id"]) %>% 
      dplyr::select(-id)
    
    # calculate correlation matrix and pvalues
    cor_mat <- my_data %>% 
      as.matrix() %>% 
      Hmisc::rcorr(type=method)
    out=flattenCorrMatrix(cor_mat$r, cor_mat$P) %>% 
      dplyr::mutate(type="cor") %>% 
      dplyr::rename("dist"="cor", "cut_p"="p") 
    return(out)
})


