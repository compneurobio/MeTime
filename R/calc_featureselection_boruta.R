#' Function to calculate dependent variables
#' @description Calculate feature importance of predictors onto a response vector using the Boruta algorithm.
#' @param object a S4 object of the class "metime_analyzer".
#' @param which_x a character defining the name of the dataset containing the predictors.
#' @param which_y a character defining the name of the dataset containing the response. Can be factor for classification or numeric vector for regression.
#' @param verbose a logical whether the progress should be printed. Default set to TRUE.
#' @param cols_for_meta_x A named list of a character vector to define columns for meta info of which_x
#' Names of the character vector in this case should always start with "color_".
#' @param cols_for_meta_y a named list of a character vector to define columns for meta info of which_y
#' @param save_per_run (Experimental) a logical on whether results should be saved in the same directory/tmp as single csv.
#' @param run_index (Experimental) a vector of runs corresponding the row number of a 
#' @param num_cores numeric input to define the number of cores that you want to use for parallel computing. Default is set to NULL which is parallel::detectCores() -1.
#' @details 
#' The Boruta algorithm is a feature selection method that compares the importance of each variable to that of randomly generated shadow variables. It uses an iterative two step approach:
#' 1. the algorithm creates shadow variables that are random permutations of the original variables. 
#' 2. Secondly, it fits a machine learning model to the original and shadow variables, and calculates the importance of each variable based on the difference between the model performance with and without the variable. 
#' 
#' Variables significantly more important than their shadow counterparts are considered important and selected for the final model. This process is repeated until all variables have been evaluated. The Boruta algorithm is particularly useful for datasets with many variables and noisy data, and can help improve the accuracy and interpretability of predictive models.
#' @seealso [Boruta::Boruta]
#' @return List of conservation index results
#' @export
#' 
setGeneric("calc_featureselection_boruta", function(object, which_x,which_y, verbose=TRUE, name="boruta_1", cols_for_meta_x=NULL, cols_for_meta_y=NULL, maxRuns=100, num_cores=NULL,save_per_run=F,run_index=NULL) standardGeneric("calc_featureselection_boruta"))
setMethod("calc_featureselection_boruta", "metime_analyser", function(object, which_x,which_y, verbose=TRUE, name="boruta_1", cols_for_meta_x=NULL, cols_for_meta_y=NULL, maxRuns=100, num_cores=NULL,save_per_run=F,run_index=NULL){
  # check variable if analysis should be run
  run=T
  
  # check if results should be saved per run - this doesn't work on server because the path cannot be extracted like this
  # This is the new approach
  if(save_per_run && !"tmp" %in% list.dirs(getwd(), recursive=F, full.names=F)) {
    dir.create(file.path(getwd(), "tmp"), showWarnings=F)
  }
  
  # sanity checks
  if(!all(length(which_x)==1, length(which_y)==1)) {
    ## check that only one character was used for which_x and which_y
    warning("calc_featureselection_boruta() can only processs two datasets at a time. If you want to use multiple datasets then merge them before using this function. Currently, exiting without making any changes")
    run=F
    out <- object
  } 
  
  if(!all(c(which_x, which_y) %in% names(object@list_of_row_data))){
    ## check that the data sets are available
    warning(paste0("calc_featureselection_boruta() could not find data set: ", setdiff(c(which_x, which_y),names(object@list_of_row_data)),"."))
    run=F
    out <- object
  }
  
  if(!is.null(cols_for_meta_x)) {
    if(grep("color_", names(cols_for_meta_x[[1]])) %>% length() != names(cols_for_meta_x[[1]])) {
      warning("calc_featureselection_boruta(): cols_for_meta_x are not named correctly. Exiting without performing the calculation")
      run=F
      out <- object
    }
  }
  
  if(!is.null(cols_for_meta_y)) {
    if(grep("color_", names(cols_for_meta_y[[1]])) %>% length() != names(cols_for_meta_y[[1]])) {
      warning("calc_featureselection_boruta(): cols_for_meta_y are not named correctly. Exiting without performing the calculation")
      run=F
      out <- object
    }
  }
  
  ## check if a dataset is numeric but should be used as factor (classification)
  my_unique_x <- apply(object@list_of_data[[which_x]], 2, function(x) length(unique(x)))%>% unique()
  my_type_x <- apply(object@list_of_data[[which_x]], 2, function(x) is.numeric(x)) %>% unique()
  if(my_unique_x <=6 && all(my_type_x == T)){
    warning(paste0("calc_featureselection_boruta() changed ", which_x, " to data frame with factors for classification"))
    object@list_of_data[[which_x]] <- object@list_of_data[[which_x]] %>% mutate_all(as.factor)
  }
  
  # run Boruta
  if(run){
    ## prepare data (join columns to be used)
    this_data <- object@list_of_data[[which_x]] %>% 
      dplyr::mutate(id =rownames(.[])) %>%
      dplyr::inner_join(object@list_of_data[[which_y]] %>% 
                          dplyr::mutate(id =rownames(.[])),
                        by="id") %>% 
      na.omit()
    this_x <- object@list_of_col_data[[which_x]]$id
    this_y <- object@list_of_col_data[[which_y]]$id
    file_path = getwd()
    out_path = file.path(getwd(), "tmp")
    
    # Set up parallel processing
    if(Sys.info()["sysname"] == "Windows"){
      max_cores <- parallel::detectCores() - 1
      if (is.null(num_cores) || num_cores>max_cores) num_cores <- max_cores # num_cores can only be max number of cores -1 
      cl <- parallel::makeCluster(num_cores)
      doParallel::registerDoParallel(cl)
      ## export data needed in analysis
      parallel::clusterExport(cl=cl, 
                              varlist=c("this_data","this_x","this_y","save_per_run","file_path", "out_path"),
                              envir = environment())
    }else{
      max_cores <- parallel::detectCores() - 1
      if (is.null(num_cores) || num_cores>max_cores) num_cores <- max_cores # num_cores can only be max number of cores -1 
      cl <- num_cores
    }
    
    
    all_runs <- 1:ncol(object@list_of_data[[which_y]])
    if(!is.null(run_index)) all_runs <- as.numeric(run_index)
    
    opb <- pbapply::pboptions(title="Running feature selection: ", type="timer")
    on.exit(pbapply::pboptions(opb))                   
    out_boruta <- pbapply::pblapply(all_runs,
                                    cl=cl,
                                    function(x){
                                      require(tidyverse)
                                      # if files are already existant
                                      if(save_per_run & paste0(this_y[x],".rds") %in% list.files(out_path)) {
                                        my_stats <- readRDS(file=file.path(out_path, paste0(this_y[x],".rds")))
                                      }
                                      else if(save_per_run & !paste0(this_y[x],".rds") %in% list.files(out_path)) {
                                        my_model <- Boruta::Boruta(x = this_data[,this_x],
                                                                   y = this_data[,this_y[x]],
                                                                   doTrace=0,
                                                                   holdHistory=T
                                        )
                                        
                                        my_stats=Boruta::attStats(my_model) %>%
                                          as.data.frame() %>%
                                          dplyr::mutate(id=this_y[x],
                                                        y=meanImp,
                                                        id_met=rownames(.[]))
                                        saveRDS(my_stats, file=file.path(out_path, paste0(this_y[x],".rds")))
                                      }else{
                                        my_model <- Boruta::Boruta(x = this_data[,this_x],
                                                                   y = this_data[,this_y[x]],
                                                                   doTrace=0,
                                                                   holdHistory=T
                                        )
                                        
                                        my_stats=Boruta::attStats(my_model) %>%
                                          as.data.frame() %>%
                                          dplyr::mutate(id=this_y[x],
                                                        y=meanImp,
                                                        id_met=rownames(.[]))
                                      }
                                      my_stats
                                    })
    out_results <- out_boruta %>% 
      plyr::rbind.fill() %>% 
      dplyr::mutate(med=id_met,met=id,
                    selected=ifelse(meanImp >= (mean(abs(meanImp),na.rm=T) + 4*sd(abs(meanImp),na.rm=T)),T,F),
                    thresh=(mean(abs(meanImp),na.rm=T) + 4*sd(abs(meanImp),na.rm=T))) %>% 
      dplyr::select(met,med,meanImp, selected,thresh)
    saveRDS(out_results, file=file.path(out_path, paste0("All_results_list.rds")))
    out_list=list()
    out_list[[paste0(name)]]<-out_results
    out <- get_make_results(object = object, 
                            data = out_list, 
                            metadata = NULL, 
                            calc_type = "feature_selection",
                            calc_info = "feature_selection",
                            name = name) %>%
      add_function_info(function_name = "calc_featureselection_boruta", 
                        params = list(which_x=which_x,
                                      which_y=which_y, 
                                      verbose=verbose, 
                                      name=name, 
                                      cols_for_meta_x=cols_for_meta_x, 
                                      cols_for_meta_y=cols_for_meta_y, 
                                      num_cores=num_cores,
                                      save_per_run=save_per_run,
                                      maxRuns=maxRuns,
                                      run_index=run_index))
  }
  
  
  if(Sys.info()["sysname"] == "Windows")  on.exit(parallel::stopCluster(cl))
  
  return(out)
})
