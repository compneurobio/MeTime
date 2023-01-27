
#' Function to calculate dependent variables
#' @description An S4 method to be applied on the metime_analyser object so as to calculate dependent variables
#' @param object An object of class metime_analyser
#' @param which_x Name of the dataset to be used for training
#' @param which_y Name of the dataset to be used for testing
#' @param output_loc path to the parent directory where in the out file wíll be stored
#' @param file_name name of the out file
#' @param verbose Information provided on steps being processed
#' @param cols_for_meta A list of a Character vector to define column names that are to be used for plotting purposes
#' @return List of conservation index results
#' @export
#' 
#' 
setGeneric("calc_featureselection_boruta", function(object, which_x, which_y, verbose, output_loc, file_name, cols_for_meta) standardGeneric("calc_featureselection_boruta"))
setMethod("calc_featureselection_boruta", "metime_analyser", function(object, which_x,which_y, verbose=F, output_loc=getwd(), file_name="boruta", cols_for_meta) {
  
  # validate arguments
  stopifnot(length(which_x)==1, length(which_y)==1, 
    all(c(which_x,which_y) %in% object@list_of_data))
  
  x_data <- object@list_of_data[[which_x]]
  y_data <- object@list_of_data[[which_y]]
  
  x_data <- x_data[intersect(rownames(x_data), rownames(y_data)),]
  y_data <- y_data[rownames(x_data),]
  
  results <- parallel::mclapply(sample(colnames(y_data)), function(i) {
    if(verbose) cat(i, "; ")
    if(paste0(file_name, ".rds") %in% list.files(output_loc)) my_results=readRDS(file=paste0(output_loc, "/", file_name, ".rds"))
      else my_results = data.frame()
      
    my_model = Boruta::Boruta(x=x_data, y=y_data[,i],
                              pValue=0.01,
                              doTrace=0,
                              mcAdj=TRUE)
    my_stats=Boruta::attStats(my_model) %>%
      as.data.frame() %>% 
      dplyr::mutate(id=i,
                    y=meanImp,
                    id_met=rownames(.[]))
   
    saveRDS(object=rbind(my_results, my_stats),file=paste0(output_loc, "/", file_name, ".rds"))
    return(my_stats) 
  }, max.cores=4) %>% do.call(what=rbind.data.frame)

  #change the below part to lapply
  count <- 1
  for(i in colnames(y_data)) {
    if(i %in% my_results$id_x) {
      dummy <- my_results$id_y[results$id_x %in% i]
      dummy <- paste(dummy, collapse="###")
      final <- data.frame()
      final$id[count] <- i
      final$covarites[count] <- dummy
      count <- count +1
    } else {
      dummy <- NA
      final$id[count] <- i
      final$covariates[count] <- dummy
      count <- count +1
    }
  }
  final <- final[order(final$id), ]
  object@list_of_col_data[[which_y]] <- object@list_of_col_data[[which_y]][order(object@list_of_data[[which_y]]$id), ]
  object@list_of_col_data[[which_y]]$covariates <- final$covariates

  # Add metadata to the results
  metadata <- get_metadata_for_columns(object = object, 
                                         which_data = i, 
                                         columns = cols_for_meta, 
                                         names = c("name", "group"), 
                                         index_of_names = "id")

  object <- object %>% get_make_results(data = list(results), 
                              metadata = metadata, 
                              calc_type = rep("feature_selection", each=length(out_sum)), 
                              calc_info = paste("feature_selection", "_", sep = ""),
                              name=file_name)
  object <- object %>% add_function_info(function_name="calc_featureselection_boruta", 
      params=list(which_x=which_x, which_y=which_y, verbose=verbose, 
          output_loc=output_loc, file_name=file_name))
  
  return(object)
})


