
#' Function to perform PARAFAC analysis 
#' @description Function for S4 object of class metime_analyser to perform PARAFAC analysis
#' @param object S4 object of class metime_analyser
#' @param which_data character vector for dataset to be used
#' @param stratifications list to define stratification for the data before calculation
#' @param nfac parameter nfac for parafac(). Numeric value to define the number of factors. Default is set to 3 
#' @param ... Additional arguments to be used for the function parafac()
#' @return An object of class PARAFAC. See multiway library for more information
#' @export
setGeneric("calc_parafac", function(object, which_data, timepoints, nfac=3, ...) standardGeneric("calc_parafac"))
setMethod("calc_parafac", "metime_analyser", function(object, which_data, timepoints, nfac=3, ...) {
        
        stopifnot(names(stratifications) %in% time)

        data_list <- get_stratified_data(object=object, which_data=which_data,
                    stratifications=stratifications)
        data <- data_list[["data"]]
        row_data <- data_list[["row_data"]]

        final_data <- lapply(stratifications$time, function(x) {
                        row <- row_data %>% filter(time==x)
                        data <- data[rownames(data) %in% rownames(row), ]
                        return(data)
                })
        names(final_data) <- stratifications$time
        samples <- lapply(final_data, function(x) {
                rownames(x) %>% gsub(pattern="_[a-z|A-Z][0-9]+", replacement="") %>%
                return()
          })
        common_samples <- Reduce(intersect, samples)
        final_data <- lapply(final_data, function(x) {
                rownames(x) <- rownames(x) %>% gsub(pattern="_[a-z|A-Z][0-9]+", replacement="")
                x <- x[rownames(x) %in% common_samples, ]
                x <- x[order(rownames(x)), ]
                return(x)
          }) 
       parafac_data <- abind(final_data, along=length(final_data))
       parafac_object <- multiway::parafac(parafac_data, nfac=nfac, ...)
       rownames(parafac_object$A) <- rownames(final_data[[1]])
       rownames(parafac_object$B) <- colnames(final_data[[1]])
       rownames(parafac_object$C) <- timepoints
       return(parafac_object)
  })


