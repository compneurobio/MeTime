#' Function to perform GAMM analysis
#' @description An S4 method to perform GAMM analysis on a particular dataset. See regression_example to use a curated dataset.
#' @param Object An S4 object of class metime_analyser                      
#' @param which_data Name of the dataset to be used
#' @param verbose Information provided on steps being processed
#' @param cols_for_meta Character vector to define column names that are to be used for plotting purposes
#' @param name character vector to define the results. Should be equal to length of which_data
#' @param stratifications List to stratify data into a subset. Usage list(name=value)
setGeneric("calc_gamm", function(object, which_data, stratifications = NULL, cols_for_meta=NULL, verbose=T,name="regression_gamm", ...) standardGeneric("calc_gamm"))
setMethod("calc_gamm", "metime_analyser", function(object, # name of analyzer object
                      which_data, # merged dataset to be used
                      stratifications = NULL,
                      cols_for_meta=NULL,
                      verbose=T,
                      name="regression_gamm_1",
                      ...) {
    out=list()
    if(grep(name, names(object@results)) %>% length() >=1) {
      warning("name of the results was previously used, using a different name")
      index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
      index <- c(0:9)[grep(index, 0:9)+1]
      name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
    }
    my_data <- get_stratified_data(which_data=which_data, 
          object=object, stratifications=stratifications)
  
  # find list of metabolites
    my_met <- object@list_of_col_data[[which_data]] %>% 
      dplyr::filter(type=="met") %>% 
      dplyr::pull(id)
    
    my_formula <- lapply(unique(my_met),function(x) object@list_of_col_data[[which_data]] %>% 
                           dplyr::filter(type=="trait") %>% 
                           dplyr::mutate(met=x)) %>% 
      do.call(what=rbind.data.frame)
    
    results=parallel::mclapply(1:100, mc.cores=parallel::detectCores(all.tests = FALSE, logical = TRUE)-1, function(x) {
                                if(verbose) cat(x, " , ") # report iteration
                                 # extract data 
                                this_data <-  my_data$data %>% 
                                   dplyr::mutate(met=get(my_formula$met[x]), 
                                                 trait=get(my_formula$id[x])) %>% 
                                   dplyr::filter(!is.na(met), !is.na(trait)) %>% 
                                   dplyr::mutate(time = as.numeric(gsub(pattern="([a-zA-Z])",replacement = "",x=time)),
                                                 subject = as.numeric(gsub(pattern="([a-zA-Z])",replacement = "",x=subject)))
                                 
                                 # extract formula information
                                formula_met <- my_formula$met[x]
                                formula_trait <- my_formula$trait[x]
                                formula_cov <- my_formula$cov[x] %>% stringr::str_split(pattern="###",simplify = T) %>% as.character()
                                formula_random <-  paste0("s(",my_formula$random[x] %>% stringr::str_split(pattern="###",simplify = T) %>% as.character(), ", bs='re')")
                                if("interaction" %in% names(my_formula)) {
                                   formula_interaction <- my_formula$interaction[x] %>% stringr::str_split(pattern="###",simplify = T) %>% as.character()
                                }
                                else {
                                  formula_interaction <- NA
                                }
                      
                                formula_extras <- my_formula$extras[x] %>% stringr::str_split(pattern="###",simplify = T) %>% as.character()
                                if(length(formula_extras) ==0) formula_extras <- NA
                                cov_smoothing <- sapply(formula_cov,function(y) ifelse(length(unique(this_data[[y]]))>3, paste0("s(",y,")"), y)) 
                                 
                                 # collapse formula
                                 formula = paste0("met ~ trait", 
                                                  ifelse(all(!is.na(cov_smoothing)), paste0("+",paste0(cov_smoothing,collapse = "+")), ""),
                                                  ifelse(all(!is.na(formula_interaction)), paste0("+",paste0(formula_interaction,collapse = "+")), ""),
                                                  ifelse(all(!is.na(formula_random)), paste0("+",paste0(formula_random,collapse = "+")), ""),
                                                  ifelse(all(!is.na(formula_extras)), paste0("+",paste0(formula_extras,collapse = "+")), "")
                                 )
                                 this_model = try(
                                   silent=T,
                                   mgcv::bam(
                                     formula=as.formula(formula), 
                                     data= this_data,
                                     method="REML"
                                   )
                                 )
                                 if(all(class(this_model)!="try-error")) {
                                   results_sum <- summary(this_model)
                                   out_this_model <- data.frame(
                                     stringsAsFactors = F,
                                     met = my_formula$met[x],
                                     trait = my_formula$id[x],
                                     level = grep(x=names(results_sum$p.coeff), pattern="trait", value=T) %>% gsub(pattern="trait",replacement=""),
                                     beta = results_sum$p.coeff[grep(x=names(results_sum$p.coeff), pattern="trait", value=T)] %>% as.numeric(),
                                     se = results_sum$se[grep(x=names(results_sum$se), pattern="trait", value=T)]%>% as.numeric(),
                                     pval = results_sum$p.pv[grep(x=names(results_sum$p.pv), pattern="trait", value=T)]%>% as.numeric(),
                                     tval = results_sum$p.t[grep(x=names(results_sum$p.t), pattern="trait", value=T)]%>% as.numeric()
                                   ) %>% 
                                     cbind(as.data.frame(t(results_sum$s.table[,1])))
                                 } else {
                                   out_this_model <- data.frame(
                                     stringsAsFactors = F,
                                     met = my_formula$met[x],
                                     trait = my_formula$id[x]#,
                                   )
                                 }
                                  return(out_this_model)
                               })
    annotated_results <- plyr::rbind.fill(results)
    
    annotated_results <- annotated_results %>% 
      dplyr::mutate(x=beta, xmin=beta-abs(se), xmax=beta+abs(se), y=met)
    
    # continue here
    out <- get_make_results(object = object, 
                            data = list(annotated_results), 
                            metadata = NULL, 
                            calc_type = "regression", 
                            calc_info = paste("gamm_calculation_for_", which_data, sep=""),
                            name = name) %>%
          add_function_info(function_name = "calc_gamm", 
                             params = list(which_data = which_data, verbose = verbose, 
                                           cols_for_meta = cols_for_meta, name = name, 
                                           stratifications = stratifications))
  
  return(out)

})

