
#' Function to scale the data
#' @description Functions for scaling
#' @examples # example to apply scaling
#' object <- mod_zscore(object, which_data="name of the dataset")
#' @param object An object of class metime_analyser
#' @param which_data Name of the dataset to be used
#' @importClassesFrom metime_analyser
#' @return An object of class metime_analyser with processed data
#' @export
setGeneric("mod_trans_zscore", function(object, which_data) standardGeneric("mod_trans_zscore"))
setMethod("mod_trans_zscore", "metime_analyser", function(object, which_data) {
  #define data to be processed
  data_position <- which(names(object@list_of_data) %in% which_data)
  data_rownames <- lapply(object@list_of_data[data_position], function(a) rownames(a))
  object@list_of_data[data_position] <- lapply(object@list_of_data[data_position], function(a) {
          a <- apply(a, 2, as.numeric)
          return(a)
    })
  object@list_of_data[data_position] = lapply(object@list_of_data[data_position], scale, center=TRUE, scale=TRUE)
  object@list_of_data[data_position] = lapply(object@list_of_data[data_position], as.data.frame)
  for(i in 1:length(data_position)) {
      rownames(object@list_of_data[[data_position[i]]]) <- data_rownames[[i]] 
  }
  out <- object
  out <- add_function_info(object=out, function_name="mod_trans_zscore", params=list(which_data=which_data))
  return(out)
})


