#' Cossectional linear models (per time point)
#' @description Fit multiple linear models (lm) onto one dataset. Covariates can be added using the covariate column of the col_data (multiple covariates can be added, separated by '###') See examples for more information.
#' @param object an S4 object of class metime_analyser                      
#' @param which_data  a character defining the name of the dataset to be used.
#' @param verbose a logical on whether to print the calculation progress. Default set to FALSE.
#' @param cols_for_meta a character vector to define column names that are to be used for plotting purposes. Default set to NULL, therby not adding columns as metadata.
#' If you want automated facet wrapping option then set your new_columns as "facet_your_name"
#' @param name a character vector to define the index within the results. Should be equal to length of which_data. Default set to regression_lm_1.
#' @param threshold a character vector to define the type of threshold for significant interactions. Default set to all availabe thresholds: c("none","nominal","li","fdr","bonferroni").
#'      allowed inputs are "li", "FDR", "bonferroni" and "nominal"(cutoff p=0.05, set as Default)
#' @param stratifications list to stratify data into a subset. Usage list(name=value). Default set to NULL, thereby not performing any type of stratification.
#' @param num_cores numeric input to define the number of cores that you want to use for parallel computing. Default is set to NULL which is parallel::detectCores() -1.
#' @param timepoint time input for cross-sectional model should be the same as the value in time column in row_data.
#' @details Add details here
#' @return a S4 object of the class metime_analyzer with analysis results appended to the result section.
#' @export
setGeneric("calc_lm", function(object, which_data, stratifications = NULL, cols_for_meta=NULL, threshold=c("none","nominal","li","fdr","bonferroni"), verbose=T,name="regression_lm_1", timepoint=NULL, num_cores=NULL) standardGeneric("calc_lm"))
setMethod("calc_lm", "metime_analyser", function(object,
                                                  which_data,
                                                  stratifications = NULL,
                                                  cols_for_meta=NULL,
                                                  threshold=c("none","nominal","li","fdr","bonferroni"),
                                                  verbose=T,
                                                  name="regression_lm_1", 
                                                  timepoint=NULL,
                                                  num_cores=NULL) {
  #sanity checks
  if(!all(c("cov", "type") %in% names(object@list_of_col_data[[which_data]]))) stop("calc_lm() requires columns with covariates (named 'cov') and type")

  while(grep(name, names(object@results)) %>% length() >= 1) {
    warning("name of the results was previously used, using a different name")
    index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
    index <- c(0:9)[grep(index, 0:9)+1]
    name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
  }
  
  #stratify data
  lm_data <- get_stratified_data(which_data=which_data, 
                                  object=object, stratifications=stratifications)
 
  

  if(is.null(timepoint)) {

      message("Timepoint has to be specified. Exiting without making any changes")
      return(object)
  }

  if(length(timepoint)>1) {
      message("Timepoint has to be of length 1. Exiting without making any changes")
      return(object)
  }


  # add time and subject to dataframe
  lm_data$data <-  lm_data$data %>% 
    dplyr::select(any_of(setdiff(names(lm_data$data), c("subject","time")))) %>% 
    dplyr::mutate(id = rownames(.[])) %>%
    dplyr::left_join(by="id", lm_data$row_data[,c("id","time","subject")]) %>%
    dplyr::filter(time %in% timepoint)
  
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
  
  # runs are defined by unique covariates versus time - not needed as only single timepoint
  # my_runs <- lapply(unique(lm_data$data$time),function(x) my_formula %>% 
  #                   dplyr::distinct(cov) %>% 
  #                   dplyr::mutate(time = x)
  #                 ) %>% 
  #  do.call(what=rbind.data.frame)

  # model calculation ----
  ## add verbose processbar 
  ## changing from mclapply to parLapply
  if(is.null(num_cores)) {
    cl <- parallel::makeCluster(spec = parallel::detectCores(all.tests = FALSE, logical = TRUE)-1, type="PSOCK")
  } else {
    cl <- parallel::makeCluster(spec = num_cores, type="PSOCK")
  }
  parallel::clusterExport(cl=cl, varlist=c("my_formula", "lm_data"), envir=environment())
  opb <- pbapply::pboptions(title="Running calc_lm(): ", type="timer")
  on.exit(pbapply::pboptions(opb))
  
  results=pbapply::pblapply(cl=cl, 1:nrow(my_formula), 
                 function(x) {
                               require(magrittr)
                               # extract data 
                               this_data <-  lm_data$data %>% 
                                 dplyr::select(any_of(setdiff(names(lm_data$data), c("subject","time")))) %>% 
                                 dplyr::mutate(id = rownames(.[])) %>%
                                 dplyr::left_join(by="id", lm_data$row_data[,c("id","time","subject")]) %>% 
                                 dplyr::mutate(met=get(my_formula$met[x]), 
                                               trait=get(my_formula$trait[x])) %>% 
                                 #dplyr::filter(!is.na(met), !is.na(trait)) %>% 
                                 dplyr::mutate(time = as.numeric(gsub(pattern="([a-zA-Z])",replacement = "",x=time)),
                                               subject = as.numeric(gsub(pattern="([a-zA-Z])",replacement = "",x=subject)))
                               
                               # extract formula information
                               formula_met <- my_formula$met[x]
                               formula_trait <- my_formula$trait[x]
                               formula_cov <- my_formula$cov[x] %>% stringr::str_split(pattern="###",simplify = T) %>% as.character() %>% .[!. %in% ""]
                               
                               formula = paste0("met ~ ",
                                                ifelse(all(!is.na(formula_trait)), "trait", ""), 
                                                ifelse(all(!is.na(formula_cov)), paste0("+",paste0(formula_cov,collapse = "+")), "")
                                        )
                               
                               this_model = try(
                                  silent=T,
                                  lm(
                                    formula=formula, 
                                    data= this_data
                                  )
                               )



                               if(all(class(this_model)!="try-error")) {
                                  results_sum <- summary(this_model)

                                  out <- data.frame(
                                      met=my_formula$met[x],
                                      trait=my_formula$trait[x],
                                      level=grep(x=rownames(results_sum$coefficients), pattern="trait", value=T) %>% gsub(pattern="trait",replacement=my_formula$trait[x]),
                                      beta = results_sum$coefficients[,1][grep(x=names(results_sum$coefficients[,1]), pattern="trait", value=T)] %>% as.numeric(),
                                      se = results_sum$coefficients[,2][grep(x=names(results_sum$coefficients[,2]), pattern="trait", value=T)]%>% as.numeric(),
                                      pval = results_sum$coefficients[,4][grep(x=names(results_sum$coefficients[,4]), pattern="trait", value=T)]%>% as.numeric(),
                                      tval = results_sum$coefficients[,3][grep(x=names(results_sum$coefficients[,3]), pattern="trait", value=T)]%>% as.numeric(),
                                      formula = formula
                                    )
                               } else {
                                  out <- data.frame(
                                      met=my_formula$met[x],
                                      trait=my_formula$trait[x],
                                      level=NA, 
                                      beta=NA,
                                      se=NA, 
                                      pval=NA,
                                      tval=NA,
                                      formula=formula
                                    )
                               }
                               
                              return(out)

                        })
  
  # combine all results and do post processing - cancelled idea - now we want the results to be separate

  annotated_results <- plyr::rbind.fill(results)
  on.exit(parallel::stopCluster(cl))
  ## modify results to for pipeline
  annotated_results <- annotated_results %>% 
    dplyr::mutate(x=beta, 
                  xmin=beta-abs(se), 
                  xmax=beta+abs(se), 
                  y=met, 
                  color= "none",
                  type=ifelse(!is.na(level), level, trait))

  #out_results <- lapply(unique(results))

  # calculate the thresholds 
  thresh_bonferroni <- 0.05/length(my_met)
  eigenvals <- cor(object@list_of_data[[which_data]][,my_met], use="pairwise.complete.obs") %>%
    eigen()
  thresh_li <- 0.05/(sum(as.numeric(eigenvals$values >= 1) + (eigenvals$values - floor(eigenvals$values))))

  annotated_results <- annotated_results %>% dplyr::mutate(li_thresh=thresh_li, bonferroni_thresh=thresh_bonferroni)
  
  # split into single results files 
  out_results <- lapply(unique(annotated_results$type), function(y){
    annotated_results[which(annotated_results$type==y),] %>% 
      dplyr::mutate(qval = p.adjust(pval, method="BH")) %>% 
      dplyr::mutate(color = ifelse(pval <= 0.05, "nominal","none")) %>% 
      dplyr::mutate(color = ifelse(qval <= 0.05, "fdr",color)) %>% 
      dplyr::mutate(color = ifelse(pval <= thresh_li, "li",color)) %>% 
      dplyr::mutate(color = ifelse(pval <= thresh_bonferroni, "bonferroni",color)) %>% 
      `rownames<-`(.[,"met"])
    })

  names(out_results) <- unique(annotated_results$type)

                            
  # add cols_for_meta
  
  if(is.null(cols_for_meta)) {
    metadata <- NULL
  } else {
    metadata <- get_metadata_for_columns(object=object, which_data=which_data, columns=cols_for_meta)
  }
                            
  # make result item
  out <- get_make_results(object = object, 
                          data = out_results, 
                          metadata = metadata, 
                          calc_type = rep("regression",each=length(out_results)),
                          calc_info = paste0("lm regression for ",names(out_results)),
                          name = name) %>%
    add_function_info(function_name = "calc_lm", 
                      params = list(which_data = which_data, 
                                    verbose = verbose, 
                                    cols_for_meta = cols_for_meta,
                                    name = name, 
                                    stratifications = stratifications,
                                    threshold=threshold))

  return(out)
})


