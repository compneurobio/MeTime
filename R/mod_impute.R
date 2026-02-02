#' Imputation of missing values
#' @description This function imputes missing values using one of the implemented imputation methods.
#' @param object a S4 object of class "metime_analyser".
#' @param which_data a character defining which dataset should be imputed.
#' @param method a character defining the method for imputation. See details for more information. Default set to 'rf'
#' @param thresh a numeric defining at which missingness level metabolites will be excluded. Default set to 0.3 (30% missingness)
#' @param path_to_diagnostics a character defining the path to where a diagnostic file is being written.
#' @details additional methods are currently not implemented such as knn, pmm, norm etc. 
#' @returns S4 object of the class "metime_analyser" with mutated col_data and row_data.
#' @export
setGeneric("mod_impute", function(object, which_data, method="rf",thresh=0.3,path_to_diagnostics=NULL) standardGeneric("mod_impute"))
setMethod("mod_impute", "metime_analyser", function(object, which_data, method="rf",thresh=0.3,path_to_diagnostics=NULL) {
  out <- object
  if(length(method)>1){
    warning("mod_impute() can only apply one method per data frame")
  }else{
    # exclusion of metabolite
    exclude_metabolites <- colSums(is.na(out@list_of_data[[which_data]]))
    exclude_thresh = nrow(out@list_of_data[[which_data]]) * thresh
    
    # extract data from analyzer object
    missing_data <- out@list_of_data[[which_data]] %>% 
      dplyr::mutate_all(as.numeric)%>% # make sure all data is numeric
      dplyr::select(-any_of(names(exclude_metabolites)[which(exclude_metabolites>exclude_thresh)])) #remove metabolites with more missingness as thresh
    
    # make sure that NA values are not defined as NaN or inf
    missing_data = as.matrix(missing_data)
    missing_data[is.nan(missing_data)] <- NA
    missing_data[is.infinite(missing_data)] <- NA
    missing_data <- as.data.frame(missing_data)
    
    if(method=="rf"){
    
      # impute data
      imputed_data <- missForest::missForest(
        xmis=missing_data,
        maxiter = 5, 
        parallelize="no",
        verbose=F,
        variablewise=TRUE
      )
      
      out@list_of_data[[which_data]] <- imputed_data[["ximp"]] %>% as.data.frame()
      if(!is.null(path_to_diagnostics)){
        saveRDS(list(
          error=imputed_data$OOBerror,
          excluded=names(exclude_metabolites)[which(exclude_metabolites>exclude_thresh)]
        ),
        file=paste0(path_to_diagnostics,"mod_impute_diagnostics.rds"))
      }
    }else if(method=="mean"){
      out@list_of_data[[which_data]]  <- missing_data %>% 
        mutate(across(everything(), ~ if_else(is.na(.), mean(., na.rm = TRUE), .)))
    }else if(method=="min"){
      out@list_of_data[[which_data]]  <- missing_data %>% 
        mutate(across(everything(), ~ if_else(is.na(.), min(., na.rm = TRUE), .)))
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