

#' Function to calculate colinearity of features
#' @description Function to calculate colinearity in a dataset
#' @param object An S4 object of class metime_analyser
#' @param which_data Dataset to check for colinearity
#' @param show_all logical. True will only filter out colinear data
#' @param name character to define the name of the result
#' @param stratifications List to stratify data into a subset. Usage list(name=values)
#' @return plotter object with data for heatmap information
#' @export
setGeneric("calc_colinearity_features", function(object, which_data, name="calc_colinearity_1", stratifications=list()) standardGeneric("calc_colinearity_features"))
setMethod("calc_colinearity_features", "metime_analyser", function(object, which_data, name="calc_colinearity_1", stratifications=list()) {
      stopifnot(which_data %in% names(object@list_of_data))
      stopifnot(length(names(object@list_of_data[[which_data]]))==length(unique(names(object@list_of_data[[which_data]]))))
      if(grep(name, names(object@results)) %>% length() >=1) {
          warning("name of the results was previously used, using a different name")
          index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
          index <- c(0:9)[grep(index, 0:9)+1]
          name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
      }
      if(!is.null(stratifications) && length(stratifications) >= 1) {
        data_list <- get_stratified_data(object=object, which_data=which_data, stratifications=stratifications)
        data <- data_list[["data"]]
      } else {
        data <- object@list_of_data[[which_data]]
      }
      numeric_cols <- names(data)[sapply(data, is.numeric) | sapply(data, is.integer)]
      if(length(numeric_cols) < 2) {
          warning("calc_colinearity_features(): need at least two numeric features for colinearity checks")
          return(object)
      }
      cor_mat <- stats::cor(data[, numeric_cols, drop=FALSE], use="pairwise.complete.obs", method="spearman")
      ut <- upper.tri(cor_mat)
      out <- data.frame(
        med_class_1=rownames(cor_mat)[row(cor_mat)[ut]],
        med_class_2=rownames(cor_mat)[col(cor_mat)[ut]],
        chi_statistic=NA_real_,
        Cramers_V=abs(cor_mat[ut]),
        cor=cor_mat[ut],
        stringsAsFactors=FALSE
      ) %>%
        dplyr::mutate(colinear=ifelse(Cramers_V > .5, TRUE, FALSE))

     
      out <- get_make_results(object=object, data=list(out), metadata=metadata, calc_type="colinearity", 
                  calc_info=paste("colinearity test of: ", which_data, sep=""), name=name)
      out <- out %>% add_function_info(function_name="calc_colinearity_features", 
                  params=list(which_data=which_data, stratifications=stratifications)) 
      return(out)
          
  })