
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
#' @param name name of the results should be of length=1
#' @param stratifications List to stratify data into a subset. Usage list(name=value)
#' @return data.frame with pairwise results
#' @export
setGeneric("calc_correlation_pairwise", function(object, which_data, method, name="calc_correlation_pairwise_1", stratifications) standardGeneric("calc_correlation_pairwise"))
setMethod("calc_correlation_pairwise", "metime_analyser", function(object, which_data, method="pearson", name="calc_correlation_pairwise_1", stratifications){
  stopifnot(all(which_data %in% names(object@list_of_data)))
  if(grep(name, names(object@results)) %>% length() >=1) {
    warning("name of the results was previously used, using a different name")
    index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
    index <- c(0:9)[grep(index, 0:9)+1]
    name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
  }
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
  if(length(stratifications)>=1) {
        dummy_data <- object@list_of_data[[which_data[1]]]
        row_data <- object@list_of_row_data[[which_data[1]]]
        stratifications <- lapply(names(stratifications), function(x) {
              row_data <- row_data[row_data[,x] %in% stratifications[[x]], ]
              return(stratifications[[x]]) 
          })
        my_data <- my_data[rownames(my_data) %in% rownames(row_data), ]
  }
    # calculate correlation matrix and pvalues
    cor_mat <- my_data %>% 
      as.matrix() %>% 
      Hmisc::rcorr(type=method)
    out=flattenCorrMatrix(cor_mat$r, cor_mat$P) %>% 
      dplyr::mutate(type="cor") %>% 
      dplyr::rename("dist"="cor", "cut_p"="p")

    out <- get_make_results(object=object, data=list(pairwise_correlation=out), metadata=NULL, calc_type="pairwise_correlation", 
                      calc_info = paste(which_data, "_and_" , method, "_pairwise_correlation", sep=""), 
                      name=name) %>%
          add_function_info(function_name="calc_correlation_pairwise",
                params=list(which_data=which_data, method=method))
    return(out)
})


