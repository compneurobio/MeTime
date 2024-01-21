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
#' @details Add details here
#' @return a S4 object of the class metime_analyzer with analysis results appended to the result section.
#' @export
setGeneric("calc_lm", function(object, which_data, stratifications = NULL, cols_for_meta=NULL, threshold=c("none","nominal","li","fdr","bonferroni"), verbose=T,name="regression_lm_1") standardGeneric("calc_lm"))
setMethod("calc_lm", "metime_analyser", function(object,
                                                  which_data,
                                                  stratifications = NULL,
                                                  cols_for_meta=NULL,
                                                  threshold=c("none","nominal","li","fdr","bonferroni"),
                                                  verbose=T,
                                                  name="regression_lm_1") {
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
  
  # add time and subject to dataframe
  lm_data$data <-  lm_data$data %>% 
    dplyr::select(any_of(setdiff(names(lm_data$data), c("subject","time")))) %>% 
    dplyr::mutate(id = rownames(.[])) %>%
    dplyr::left_join(by="id", lm_data$row_data[,c("id","time","subject")])

  if(length(unique(lm_data$data$time)) > 1) {
      message("There are more than one timepoints. Exiting without making any changes")
      return(object)
  }
  
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
  #my_runs <- lapply(unique(lm_data$data$time),function(x) my_formula %>% 
  #                   dplyr::distinct(cov) %>% 
  #                   dplyr::mutate(time = x)
  #                 ) %>% 
  #  do.call(what=rbind.data.frame)

  
  
  results=lapply(1:nrow(my_formula), 
                 function(x) {
                               # extract data 
                               met = unique(my_formula$met[which(my_formula$cov==my_runs$cov[x])])
                               trait = unique(my_formula$trait[which(my_formula$cov==my_runs$cov[x])])
                               cov = my_runs$cov[x] %>% stringr::str_split(pattern="###",simplify = T) %>% as.character()
                               cov = setdiff(cov, c("time","subject")) # time cannot be a covariate
                               
                               # get full data
                               this_data <- lm_data$data %>% 
                                 dplyr::filter(time==my_runs$time[x]) %>% 
                                 dplyr::select(any_of(met), any_of(trait),any_of(cov)) #%>% 
                                 #dplyr::mutate(across(trait, as.numeric)) %>% 
                                 #dplyr::filter(!is.na(trait), !is.na(met))
                               
                              met_data <- MatrixEQTL::SlicedData$new()
                              met_data$CreateFromMatrix(t(this_data %>% 
                                                            dplyr::select(any_of(met))%>% 
                                                            dplyr::mutate_all(as.numeric))) # leave in features X metabolites
                                 
                              t_data <- MatrixEQTL::SlicedData$new()
                              t_data$CreateFromMatrix(t(this_data %>%  
                                                          dplyr::select(any_of(trait)) %>% 
                                                          dplyr::mutate_all(as.numeric)
                                                           ))
                                 cov_data <- MatrixEQTL::SlicedData$new()
                                 cov_data$CreateFromMatrix(t(this_data %>% 
                                                               dplyr::select(any_of(cov))
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
                                     dplyr::rename(trait=gene, met=snps, pval=pvalue) %>%
                                     dplyr::mutate(x=beta, y=met)
                                   out$time = as.character(my_runs$time[x])
                                   out$type = out$time
                                 }else{
                                   out <- NULL
                                 }
                                 out
                             })
  
  # combine all results and do post processing - cancelled idea - now we want the results to be separate

  out_results <- results %>%
    plyr::rbind.fill()

  #out_results <- lapply(unique(results))
  
  all_results <- lapply(unique(out_results$trait), function(a) {
      this_results <- out_results %>% dplyr::filter(trait %in% a)
      eigenvals <- cor(lm_data$data[ ,my_met],use = 'pairwise.complete.obs') %>% eigen()
      li_thresh <- 0.05/(sum(as.numeric(eigenvals$values >= 1) + (eigenvals$values - floor(eigenvals$values))))

      res <- this_results %>%
        dplyr::mutate(qval = p.adjust(pval, method="BH")) %>% 
        dplyr::mutate(color= "none") %>% ## add color layers
        dplyr::mutate(color=ifelse(pval<=0.05, "nominal",color))%>% 
        dplyr::mutate(color=ifelse(pval<=li_thresh, "li",color)) %>%
        dplyr::mutate(li_thresh=li_thresh) %>% 
        dplyr::mutate(color=ifelse(pval<=0.05/length(my_met), "bonferroni",color)) %>%
        dplyr::mutate(bonferroni_thresh=0.05/length(my_met))
      return(res)
  })

  names(all_results) <- unique(out_results$trait)




  # for(x_trait in unique(out_results$trait)){
  #   this_results <- out_results[which(out_results$time == x_time),]
  #   eigenvals <- cor(lm_data$data[which(lm_data$data$time == x_time),my_met],use = 'pairwise.complete.obs') %>%eigen()
  #   li_thresh <- 0.05/(sum(as.numeric(eigenvals$values >= 1) + (eigenvals$values - floor(eigenvals$values))))
    
  #   all_results[[paste0("t",x_time)]] <- lapply(unique(this_results$trait), function(x) {
  #     res <- this_results[which(this_results$trait == x),] %>% 
  #       dplyr::mutate(qval = p.adjust(pval, method="BH")) %>% 
  #       dplyr::mutate(color= "none") %>% ## add color layers
  #       dplyr::mutate(color=ifelse(pval<=0.05, "nominal",color))%>% 
  #       dplyr::mutate(color=ifelse(pval<=li_thresh, "li",color)) %>%
  #       dplyr::mutate(li_thresh=li_thresh) %>% 
  #       dplyr::mutate(color=ifelse(pval<=0.05/length(my_met), "bonferroni",color)) %>%
  #       dplyr::mutate(bonferroni_thresh=0.05/length(my_met))  
  #   }) %>% 
  #     plyr::rbind.fill()
  # }

  # fixing the bug here
  #all_results <- lapply(seq_along(all_results), function(x) {
  #    rownames(all_results[[x]]) <- all_results[[x]]$met
  #    return(all_results[[x]])
  #})
                            
  # add cols_for_meta
  
  if(is.null(cols_for_meta)) {
    metadata <- NULL
  } else {
    metadata <- get_metadata_for_columns(object=object, which_data=which_data, columns=cols_for_meta)
  }
                            
  # make result item
  out <- get_make_results(object = object, 
                          data = all_results, 
                          metadata = metadata, 
                          calc_type = rep("regression",each=length(all_results)),
                          calc_info = paste0("lm regression for ",names(all_results)),
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
