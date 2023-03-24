#' Calculate significance threshold after correcting for independent tests. 
#' @description Eigenvalue based calculation of independent test (removing colinearity) described by Li et al. (See details for more)
#' @param object a S4 object of the class "metime_analyzer".
#' @param which_data a character to define which dataset is to be used.
#' @param verbose a logical to print out the number of independent tests. Default is set to FALSE.
#' @details A detailed description of the method can be found in \href{https://doi.org/10.1007%2Fs00439-011-1118-2}{Li et al. 2012, Evaluating the effective numbers of independent tests and significant p-value thresholds in commercial genotyping arrays and public imputation reference datasets}
#' @return a numeric value describing the Li threshold.
#' @export
setGeneric("get_li_thresh", function(object, which_data, verbose=F) standardGeneric("get_li_thresh"))
setMethod("get_li_thresh", "metime_analyser", function(object, which_data, verbose=F) {
        data <- object@list_of_data[[which_data]]
        data <- data %>% as.matrix() %>% .[,] %>% as.data.frame()  
        cordat <- cor(data)
        eigenvals <- eigen(cordat)$values
        li.thresh <- sum(as.numeric(eigenvals >= 1) + (eigenvals - floor(eigenvals)))
        if(verbose) {
          cat(paste("The number of independent tests are: ", li.thresh))
          cat("\n")
        }
        out <- 0.05/li.thresh
        return(out)   
  }) 

