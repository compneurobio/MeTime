
#' Function to calculate dependent variables
#' @description An S4 method to be applied on the metime_analyser object so as to calculate dependent variables
#' @param object An object of class metime_analyser
#' @param which_x Name of the dataset to be used for training
#' @param which_y Name of the dataset to be used for testing
#' @param output_loc path to the parent directory where in the out file wíll be stored
#' @param file_name name of the out file
#' @param verbose Information provided on steps being processed
#' @param cols_for_meta_x A list of a Character vector to define columns for meta info of which_x
#' Names of the character vector in this case should always start with "color_".
#' @param cols_for_meta_y A list of a Character vector to define columns for meta info of which_y
#' @return List of conservation index results
#' @export
#' 
#' 
setGeneric("calc_featureselection_boruta", function(object, which_x, which_y, verbose=TRUE, output_loc, file_name, cols_for_meta_x=NULL, cols_for_meta_y=NULL) standardGeneric("calc_featureselection_boruta"))
setMethod("calc_featureselection_boruta", "metime_analyser", function(object, which_x,which_y, verbose=TRUE, output_loc=getwd(), file_name="boruta", cols_for_meta_x=NULL, cols_for_meta_y=NULL) {
  
  # validate arguments and sanity checks
  if(!all(length(which_x)==1, length(which_y)==1)) {
    warning("Only two datasets can be used at a time. If you want to use multiple datasets then merge them before using this function. Currently, exiting without making any changes")
    return(object)
  } 
  if(!all(c(which_x,which_y) %in% names(object@list_of_data))) {
    warning("The datasets mentioned are not found in the metime_analyser object. Exiting without making any changes.")
    return(object)
  }
  if(!is.null(cols_for_meta_x)) {
    if(grep("color_", names(cols_for_meta_x[[1]])) %>% length() != names(cols_for_meta_x[[1]])) {
      warning("The cols_for_meta_x are not named correctly. Exiting without performing the calculation")
      return(object)
    }
  } else {
    warning("cols_for_meta_x is empty. Continuing with generating results")
  }
  
  # extracting data of interest and making sure they have the same samples
  x_data <- object@list_of_data[[which_x]]
  y_data <- object@list_of_data[[which_y]]
  x_data <- x_data[intersect(rownames(x_data), rownames(y_data)),]
  y_data <- y_data[rownames(x_data),]
  
  #Using parLappy because mclapply doesn't work for windows
  cl <- parallel::makeCluster(parallel::detectCores(all.tests = FALSE, logical = TRUE)-1)
  parallel::clusterExport(cl=cl, varlist=c("file_name", "output_loc", "x_data", "y_data"), envir=environment())
  opb <- pbapply::pboptions(title="Running calc_featureselection_boruta(): ", type="timer")
  on.exit(pbapply::pboptions(opb))
  on.exit(parallel::stopCluster(cl))
  results <- pbapply::pblapply(cl=cl, sample(colnames(y_data)), function(i) {
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
   
    saveRDS(object=rbind.data.frame(my_results, my_stats),file=paste0(output_loc, "/", file_name, ".rds"))
    return(my_stats) 
  }) %>% do.call(what=rbind.data.frame)

  final <- lapply(colnames(y_data), function(a) {
        if(a %in% my_results$id_y) {
            dummy <- my_results$id_x[results$id_y %in% a]
            dummy <- paste(dummy, collapse="###")
            return(data.frame(id=a, covariates=dummy))
        } else {
            return(data.frame(id=a, covariates=NA))
        }
    }) %>% do.call(what=rbind.data.frame)
  final <- final[order(final$id), ]
  object@list_of_col_data[[which_y]] <- object@list_of_col_data[[which_y]][order(object@list_of_data[[which_y]]$id), ]
  object@list_of_col_data[[which_y]]$covariates <- final$covariates

  # Add metadata to the results
  metadata_y <- get_metadata_for_columns(object = object, 
                                         which_data = which_y, 
                                         columns = cols_for_meta_y 
                                         )
  metadata_x <- get_metadata_for_columns(object=object,
                                          which_data=which_x,
                                          columns=cols_for_meta_x
                                          )

  object <- object %>% get_make_results(data = list(feature_selection=results), 
                              metadata = NULL, 
                              calc_type = "feature_selection", 
                              calc_info = paste("feature_selection", "_", sep = ""),
                              name=file_name)
  object <- object %>% add_function_info(function_name="calc_featureselection_boruta", 
      params=list(which_x=which_x, which_y=which_y, verbose=verbose, 
          output_loc=output_loc, file_name=file_name))

  results <- object@results[[file_name]]

  data <- results$plot_data[[1]]

  colnames(metadata_x)[colnames(metadata_x) %in% "id"] <- "id_x"
  colnames(metadata_y)[colnames(metadata_y) %in% "id"] <- "id_y"  
  data <- data %>%
            dplyr::left_join(metadata_y, by="id_y")
  data <- data %>% 
            dplyr::left_join(metadata_x, by="id_x")

  results$plot_data[[1]] <- data
  object@results[[file_name]] <- results
  return(object)
})


