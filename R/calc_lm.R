#' Calculation of crossectional linear models (per time point)
#' @description Fits multiple linear models (LM).
#' @param object an S4 object of class metime_analyser                      
#' @param which_data  a character defining the name of the dataset to be used.
#' @param verbose a logical on whether to print the calculation progress. Default set to FALSE.
#' @param cols_for_meta a character vector to define column names that are to be used for plotting purposes. Default set to NULL, therby not adding columns as metadata.
#' @param name a character vector to define the index within the results. Should be equal to length of which_data. Default set to regression_lm_1.
#' @param threshold a character of length 1 to define the type of threshold for significant interactions. 
#'      allowed inputs are "li", "FDR", "bonferroni" and "nominal"(cutoff p=0.05, set as Default)
#' @param stratifications list to stratify data into a subset. Usage list(name=value). Default set to NULL, thereby not performing any type of stratification.
#' @details Add details here
#' @return a S4 object of the class metime_analyzer with analysis results appended to the result section.
#' @export
setGeneric("calc_lm", function(object, which_data, stratifications = NULL, cols_for_meta=NULL, threshold="nominal", verbose=T,name="regression_lm_1") standardGeneric("calc_lm"))
setMethod("calc_lm", "metime_analyser", function(object,
                                                  which_data,
                                                  stratifications = NULL,
                                                  cols_for_meta=NULL,
                                                  threshold="nominal",
                                                  verbose=T,
                                                  name="regression_lm_1") {
  out=list()
  if(grep(name, names(object@results)) %>% length() >= 1) {
    warning("name of the results was previously used, using a different name")
    index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
    index <- c(0:9)[grep(index, 0:9)+1]
    name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
  }
  
  #stratify data
  lm_data <- get_stratified_data(which_data=which_data, 
                                  object=object, stratifications=stratifications)
  stopifnot(all(c("cov", "type") %in% names(lm_data$col_data)))
  
  # add time and subject to dataframe
  lm_data$data <-  lm_data$data %>% 
    dplyr::select(any_of(setdiff(names(lm_data$data), c("subject","time")))) %>% 
    dplyr::mutate(id = rownames(.[])) %>%
    dplyr::left_join(by="id", lm_data$row_data[,c("id","time","subject")])
  
  # find list of metabolites
  my_met <- object@list_of_col_data[[which_data]] %>% 
    dplyr::filter(type=="met") %>% 
    dplyr::pull(id)
  
  # find list of traits
  my_trait <- object@list_of_col_data[[which_data]] %>% 
    dplyr::filter(type=="trait") %>% 
    dplyr::pull(id)

  # find all formulas
  my_formula <- lapply(unique(my_trait),function(x)
    test <- object@list_of_col_data[[which_data]] %>% 
      dplyr::select(cov, type, id) %>% 
      dplyr::filter(id %in% my_met) %>% 
      dplyr::rename("met"="id") %>% 
      dplyr::mutate(trait=x,
                    cov = paste0(ifelse(is.na(cov), "", cov),
                                 object@list_of_col_data[[which_data]]$cov[which(object@list_of_col_data[[which_data]]$id==x)]))) %>% 
    do.call(what=rbind.data.frame)
  
  # prep for parallel computation
  my_runs <-lapply(unique(lm_data$data$time),function(x) my_formula %>% 
                     dplyr::distinct(cov) %>% 
                     dplyr::mutate(time = x)
                   ) %>% 
    do.call(what=rbind.data.frame)
  
  if(verbose){
    add_progress_bar <- txtProgressBar(min = 1, max = nrow(my_runs), style = 3)
    on.exit(close(add_progress_bar))
    cl <- parallel::makeCluster(parallel::detectCores(all.tests = FALSE, logical = TRUE)-1)
    on.exit(parallel::stopCluster(cl))
  }
  
  results=parallel::mclapply(1:nrow(my_runs), 
                             mc.cores=parallel::detectCores(all.tests = FALSE, logical = TRUE)-1,
                             mc.preschedule = TRUE,
                             function(x) {
                              # if(verbose) cat(x, " , ") # report iteration
                               # extract data 
                               met = unique(my_formula$met[which(my_formula$cov==my_runs$cov[x])])
                               trait = unique(my_formula$trait[which(my_formula$cov==my_runs$cov[x])])
                               cov = my_runs$cov[x] %>% stringr::str_split(pattern="###",simplify = T) %>% as.character()
                               cov = setdiff(cov, c("time","subject")) # time cannot be a covariate
                               
                               this_data <- lm_data$data %>% 
                                 dplyr::filter(time==my_runs$time[x]) %>% 
                                 dplyr::select(all_of(c(met, trait,cov))
                                 )
                               for(i in trait){
                                 this_data <- this_data %>% 
                                   dplyr::filter(!is.na(get(i)))
                                 
                                 this_data[,i] = as.numeric(this_data[,i])
                               }
                               
                              met_data <- MatrixEQTL::SlicedData$new()
                              met_data$CreateFromMatrix(t(this_data %>% dplyr::select(all_of(met)))) # leave in features X metabolites
                                 
                              t_data <- MatrixEQTL::SlicedData$new()
                              t_data$CreateFromMatrix(t(this_data %>%  dplyr::select(all_of(trait))
                                                           ))
                                 cov_data <- MatrixEQTL::SlicedData$new()
                                 cov_data$CreateFromMatrix(t(this_data %>% 
                                                               dplyr::select(all_of(cov))
                                                             ))
                                 
                                 this_model = try(
                                   MatrixEQTL::Matrix_eQTL_main(
                                   snps = met_data, #met_data
                                   gene = t_data, #t_data
                                   cvrt = cov_data,
                                   pvOutputThreshold = 1, # output all p-values
                                   verbose = F
                                   ),
                                 silent = T)
                                 if(all(class(this_model)!="try-error")){
                                   out <- this_model$all$eqtls %>%
                                     dplyr::rename(trait=gene, met=snps) %>%
                                     dplyr::mutate(x=beta, y=met)
                                   out$time = as.character(my_runs$time[x])
                                   out$type = out$time
                                 }else{
                                   out <- my_formula %>% 
                                     dplyr::select(met,trait) %>% 
                                     dplyr::mutate(statistic=NA,pvalue=NA,FDR=NA,beta=NA, x=NA,y=met)
                                   out$time = as.character(my_runs$time[x])
                                   out$type = out$time
                                 }
                                 rownames(out) <- out$y
                                 out
                               })
  if(verbose) setTxtProgressBar(add_progress_bar, nrow(my_runs))
  names(results) <- paste0(name,"_" ,unique(lm_data$data$time))
  
  if(is.null(cols_for_meta)) {
    metadata <- NULL
  } else {
    metadata <- get_metadata_for_columns(object=object, which_data=which_data, columns=cols_for_meta)
  }
                            
  # continue here
  out <- get_make_results(object = object, 
                          data = results, 
                          metadata = metadata, 
                          calc_type = rep("regression",each=length(results)),
                          calc_info = paste0("lm_calculation_for_", which_data,"at time ",unique(lm_data$data$time)),
                          name = name) %>%
    add_function_info(function_name = "calc_lm", 
                      params = list(which_data = which_data, 
                                    verbose = verbose, 
                                    cols_for_meta = cols_for_meta,
                                    name = name, 
                                    stratifications = stratifications))
  
  return(out)
})

