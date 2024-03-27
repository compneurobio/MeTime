#' Function to calculate fold change
#' @description Function to calculate relative log transformed concentration of 
#' metabolites with respect to baseline concentration (foldchange) of multiple
#' datasets at once
#' @param object an S4 object of class metime_analyser.                     
#' @param which_data a character vector defining the name of the dataset/s to be used. 
#' @param log2 logical input. TRUE denotes that the data is already log transformed.
#' @param col character vector to select the metabolites. Set to "all" for all metabolites.
#' @return a S4 object of the class metime_analyzer with datasets updated with foldchange values.
#' @export
#'

setGeneric("mod_trans_foldchange", function(object, which_data, log2=T, col="all") standardGeneric("mod_trans_foldchange"))

setMethod("mod_trans_foldchange", "metime_analyser", function(object, which_data, log2=T, col="all") {
    #log transform if necessary
    if(!log2){
      object <- object %>% 
        mod_trans_log(which_data=which_data, base=2)
    }
    for(ii in data_position) {
        #identify columns to transform
        if(all(col=="all")) my_cols <- object@list_of_col_data[[ii]]$id 
        else my_cols <- intersect(names(object@list_of_data[[ii]]),col)
        
        #get data ordering
        base_time <- object@list_of_row_data[[ii]]$time %>% as.numeric() %>% min(na.rm=T)
        base_id <- object@list_of_row_data[[ii]][which(object@list_of_row_data[[ii]]$time==base_time),c("id","subject","time")]
        rownames(base_id) <- base_id$subject
        # sort time points
        timepoints <-object@list_of_row_data[[ii]]$time %>% as.numeric() %>% unique() %>% sort()
        
        # calculate per time point
        new_data <- lapply(timepoints,function(yy){
          this_nonbase_id <- object@list_of_row_data[[ii]][which(object@list_of_row_data[[ii]]$time==yy),c("id","subject","time")]
          rownames(this_nonbase_id) <- this_nonbase_id$subject
          this_nonbase_id <- this_nonbase_id[intersect(rownames(this_nonbase_id),rownames(base_id)),]
          this_base_id <- base_id[rownames(this_nonbase_id),]
          out <- object@list_of_data[[ii]][this_nonbase_id$id,my_cols] - object@list_of_data[[ii]][this_base_id$id,my_cols]
          out$id <- this_nonbase_id$id
          out
        }) %>% 
          plyr::rbind.fill()
        rownames(new_data) <- new_data$id
        new_data$id <- NULL
        object@list_of_row_data[[ii]] <- object@list_of_row_data[[ii]] %>% 
            dplyr::filter(id %in% rownames(new_data))
        object@list_of_data[[ii]] <- object@list_of_data[[ii]][object@list_of_row_data[[ii]]$id,]
    
        object@list_of_data[[ii]][,my_cols] <- new_data[object@list_of_row_data[[ii]]$id,my_cols]
    }
    out <- object %>% 
    add_function_info(function_name = "mot_trans_foldchange", 
                      params = list(which_data = which_data, log2 = log2))
    return(out)

})


