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

  if(grep(name, names(object@results)) %>% length() >= 1) {
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
  
  results=parallel::mclapply(1:nrow(my_runs), 
                             mc.cores=parallel::detectCores(all.tests = FALSE, logical = TRUE)-1,
                             mc.preschedule = TRUE,
                             function(x) {
                               if(verbose) cat(x, " , ") # report iteration
                               # extract data 
                               met = unique(my_formula$met[which(my_formula$cov==my_runs$cov[x])])
                               trait = unique(my_formula$trait[which(my_formula$cov==my_runs$cov[x])])
                               cov = my_runs$cov[x] %>% stringr::str_split(pattern="###",simplify = T) %>% as.character()
                               cov = setdiff(cov, c("time","subject")) # time cannot be a covariate
                               
                               # get full data
                               this_data <- lm_data$data %>% 
                                 dplyr::filter(time==my_runs$time[x]) %>% 
                                 dplyr::select(all_of(c(met, trait,cov))) %>% 
                                 dplyr::mutate_at(trait, as.numeric) %>% 
                                 dplyr::filter(!is.na(trait), !is.na(met))
                               
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
                                     dplyr::rename(trait=gene, met=snps, pval=pvalue) %>%
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
                                 out <- out %>% dplyr::mutate(
                                   color= "none")
  
                                      for(i in intersect(c("none","nominal","li","fdr","bonferroni"),threshold)){
                                        if(i=="none"){
                                          out<-out %>% 
                                            dplyr::mutate(color="none") 
                                        }else if(i=="nominal"){
                                          out<-out %>% 
                                            dplyr::mutate(color=ifelse(pval<=0.05, "nominal",color))
                                        }else if(i == "li"){
                                          eigenvals <- cor(lm_data$data[which(lm_data$data$time == my_runs$time[x]),my_met]) %>%
                                            eigen()
                                          li_thresh <- 0.05/(sum(as.numeric(eigenvals$values >= 1) + (eigenvals$values - floor(eigenvals$values))))
                                          out<-out%>% 
                                            dplyr::mutate(color=ifelse(pval<=li_thresh, "li",color))
                                        }else if(i=="fdr"){
                                          out<-out %>% 
                                            dplyr::mutate(fdr=p.adjust(pval, method="BH")) %>% 
                                            dplyr::mutate(color=ifelse(fdr<=0.05, "fdr",color)) %>% 
                                            dplyr::select(-fdr)
                                        }else if(i=="bonferroni"){
                                          out=out %>% 
                                            dplyr::mutate(color=ifelse(pval<=0.05/length(my_met), "bonferroni",color))
                                        }
                                      }
                                 out
                               })
  
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

