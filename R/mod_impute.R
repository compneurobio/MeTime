#' Imputation of missing values
#' @description This function imputes missing values using one of the implemented imputation methods.
#' @param object a S4 object of class "metime_analyser".
#' @param which_data a character defining which dataset should be imputed.
#' @param method a character defining the method for imputation. See details for more information. Default set to 'rf'
#' @details TBD 
#' @returns S4 object of the class "metime_analyser" with mutated col_data and row_data.
#' @export
setGeneric("mod_impute", function(object, which_data, method="rf") standardGeneric("mod_impute"))
setMethod("mod_impute", "metime_analyser", function(object, which_data, method="rf") {
  out <- object
  if(length(method)>1){
    warning("mod_impute() can only apply one method per data frame")
  }else{
    if(method=="rf"){
      imputed_data <- missForest::missForest(
        xmis=out@list_of_data[[which_data]],
        maxiter = 5, 
        parallelize="no",
        verbose=F,
        variablewise=TRUE
      )
      out@list_of_data[[which_data]] <- imputed_data[["ximp"]]
    }else if(method=="mean"){
      # add
    }else if(method=="min"){
      # add
    }else if(method=="knn"){
      # add
    }else if(method=="pmm"){
      # add
    }else if(method=="norm"){
      # add
    }else if(method=="bpca"){
      # add
    }else if(method=="ppca"){
      # add
    }else if(method=="2l.pmm"){
      # add
    }else if(method=="2l.norm"){
      # add
    }
  }
  
  out <- out %>% 
		  add_function_info(function_name="mod_impute", params=list(which_data=which_data, method=method))
		return(out)
	})