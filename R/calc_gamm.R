#' Calculation of generalized additive mixed models (GAMMs)
#' @description Fits multiple generalized additive mixed models (GAMMs) to a longitudinal dataset. 
#' @param object an S4 object of class metime_analyser                      
#' @param which_data  a character defining the name of the dataset to be used.
#' @param verbose a logical on whether to print the calculation progress. Default set to FALSE.
#' @param cols_for_meta a character vector to define column names that are to be used for plotting purposes. Default set to NULL, therby not adding columns as metadata.
#' @param name a character vector to define the index within the results. Should be equal to length of which_data. Default set to regression_gamm_1.
#' @param stratifications list to stratify data into a subset. Usage list(name=value). Default set to NULL, thereby not performing any type of stratification.
#' @param random a character vector defining which variables should be treated as random effects. Default set to "subject".
#' @param interaction a character vector defining which interaction terms should be added to the model. Default set to NULL, with no interaction added.
#' @details The calculation function fits multiple generalized additive mixed models (GAMMs) on a longitudinal dataset. Here, one model fits one metabolite vs one trait. The degree of smoothness of a model term is estimated as part of the fitting. 
#' @return a S4 object of the class metime_analyzer with analysis results appended to the result section.
#' @export
setGeneric("calc_gamm", function(object, which_data, stratifications = NULL, cols_for_meta=NULL, random="subject",interaction = NULL,verbose=T,name="regression_gamm_1") standardGeneric("calc_gamm"))
setMethod("calc_gamm", "metime_analyser", function(object,
                                                   which_data,
                                                   stratifications = NULL,
                                                   cols_for_meta=NULL,
                                                   random="subject",
                                                   interaction=NULL,
                                                   verbose=T,
                                                   name="regression_gamm_1") {
  out=list()
  if(grep(name, names(object@results)) %>% length() >= 1) {
    warning("name of the results was previously used, using a different name")
    index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
    index <- c(0:9)[grep(index, 0:9)+1]
    name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
  }
  
  #stratify data
  gamm_data <- get_stratified_data(which_data=which_data, 
                                 object=object, stratifications=stratifications)
  stopifnot(all(c("cov","type") %in% names(gamm_data$col_data)))
  
  # find list of metabolites
  my_met <- object@list_of_col_data[[which_data]] %>% 
    dplyr::filter(type=="met") %>% 
    dplyr::pull(id)
  
  # find list of traits
  my_trait <- object@list_of_col_data[[which_data]] %>% 
    dplyr::filter(type=="trait") %>% 
    dplyr::pull(id)
  
  my_formula <- lapply(unique(my_trait),function(x)
      test <- object@list_of_col_data[[which_data]] %>% 
        dplyr::select(cov, type, id) %>% 
        dplyr::filter(id %in% my_met) %>% 
        dplyr::rename("met"="id") %>% 
        dplyr::mutate(trait=x,
                      cov = paste0(ifelse(is.na(cov), "", cov),
                                   object@list_of_col_data[[which_data]]$cov[which(object@list_of_col_data[[which_data]]$id==x)]))) %>% 
    do.call(what=rbind.data.frame)
  
  results=parallel::mclapply(1:nrow(my_formula), 
                             mc.cores=parallel::detectCores(all.tests = FALSE, logical = TRUE)-1, function(x) {
    if(verbose) cat(x, " , ") # report iteration
    # extract data 
    this_data <-  gamm_data$data %>% 
      dplyr::select(any_of(setdiff(names(gamm_data$data), c("subject","time")))) %>% 
      dplyr::mutate(id = rownames(.[])) %>%
      dplyr::left_join(by="id", gamm_data$row_data[,c("id","time","subject")]) %>% 
      dplyr::mutate(met=get(my_formula$met[x]), 
                    trait=get(my_formula$trait[x])) %>% 
      dplyr::filter(!is.na(met), !is.na(trait)) %>% 
      dplyr::mutate(time = as.numeric(gsub(pattern="([a-zA-Z])",replacement = "",x=time)),
                    subject = as.numeric(gsub(pattern="([a-zA-Z])",replacement = "",x=subject)))
    
    # extract formula information
    formula_met <- my_formula$met[x]
    formula_trait <- my_formula$trait[x]
    formula_cov <- my_formula$cov[x] %>% stringr::str_split(pattern="###",simplify = T) %>% as.character()
    smoothing_cov <- sapply(formula_cov,function(y) ifelse(length(unique(this_data[[y]]))>3, paste0("s(",y,")"), y)) 
    
    formula_random <-  ifelse(!is.null(random), paste0("s(",random, ", bs='re')"), NA)
    formula_interaction <- ifelse(!is.null(interaction), paste0("s(",interaction, ")"), NA)
    
    
    # collapse formula
    formula = paste0("met ~ trait", 
                     ifelse(all(!is.na(smoothing_cov)), paste0("+",paste0(smoothing_cov,collapse = "+")), ""),
                     ifelse(all(!is.na(formula_interaction)), paste0("+",paste0(formula_interaction,collapse = "+")), ""),
                     ifelse(all(!is.na(formula_random)), paste0("+",paste0(formula_random,collapse = "+")), "")
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
        trait = my_formula$trait[x],
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
        trait = my_formula$trait[x]#,
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
                          calc_type = "GAMM", 
                          calc_info = paste("gamm_calculation_for_", which_data, sep=""),
                          name = name) %>%
    add_function_info(function_name = "calc_gamm", 
                      params = list(which_data = which_data, 
                                    verbose = verbose, 
                                    cols_for_meta = cols_for_meta, 
                                    random=random,
                                    interaction=interaction,
                                    name = name, 
                                    stratifications = stratifications))
  
  return(out)
  
})

