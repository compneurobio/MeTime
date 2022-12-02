
#' Function to perform PARAFAC analysis 
#' @description Method to be applied on S4 object of class metime_analyser to perform PARAFAC analysis
#' @param object S4 object of class metime_analyser
#' @param which_data character vector for dataset to be used
#' @param timepoints character vector to define timepoints of interest
#' @param nfac parameter nfac for parafac(). Numeric value to define the number of factors. Default is set to 3 
#' @param ... Additional arguments to be used for the function parafac()
#' @return An object of class PARAFAC. See multiway library for more information
#' @export
setGeneric("calc_parafac", function(object, which_data, timepoints, nfac=3, ...) standardGeneric("calc_parafac"))
setMethod("calc_parafac", "metime_analyser", function(object, which_data, timepoints, nfac=3, ...) {
        object <- mod_split_acc_to_time(object)
        list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
        list_of_data <- lapply(list_of_data, function(x) {
              x <- x[names(x) %in% timepoints]
              x <- lapply(x, function(a) return(a[order(rownames(a)), ]))
              return(x)
          })
        final_data <- list()
        list_of_data <- unname(list_of_data)
        for(i in 1:length(timepoints)) {
            timepoint_list <- lapply(list_of_data, function(x) return(x[[i]]))
            final_data[[i]] <- do.call(cbind, timepoint_list)
        }
        names(final_data) <- timepoints
        samples <- lapply(final_data, function(x) {
                x <- rownames(x)
                x <- unlist(lapply(strsplit(x, split="_"), function(y) return(y[1])))
                return(x)
          })
        common_samples <- Reduce(intersect, samples)
        final_data <- lapply(final_data, function(x) {
                rownames(x) <- unlist(lapply(strsplit(rownames(x), split="_"), function(y) return(y[1])))
                x <- x[rownames(x) %in% common_samples, ]
                x <- x[order(rownames(x)), ]
                return(x)
          }) 
       parafac_data <- abind::abind(final_data, along=3)
       parafac_object <- multiway::parafac(parafac_data, nfac=nfac)
       rownames(parafac_object$A) <- paste(rownames(final_data[[1]]), timepoints[1], sep="_")
       rownames(parafac_object$B) <- colnames(final_data[[1]])
       rownames(parafac_object$C) <- timepoints
       return(parafac_object)
  })


