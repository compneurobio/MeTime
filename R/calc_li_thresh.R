
#' Function to calculate multiple tests using method described by li
#' @description li test to check for colinearlity and use it for feature selection
#' @param object an S4 object of class metime_analyser
#' @param which_data dataset to be used for testing
#' @param verbose Logical to print out the number of independent tests
#' @return li threshold value
#' @export
setGeneric("calc_li_thresh", function(object, which_data, verbose) standardGeneric("calc_li_thresh"))
setMethod("calc_li_thresh", "metime_analyser", function(object, which_data, verbose) {
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

