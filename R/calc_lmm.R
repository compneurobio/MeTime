#' Function to perform Linear Mixed Models
#' @description Function to perform linear mixed models on dataset of interest
#' @param object An S4 object of class metime_analyser
#' @param which_data Dataset to be used 
#' @param formula Formulae to define the equation to be used in linear mixed models
#' @param cores Numeric to define the number of cores to be used 
#' @return Returns a plotter object for plotting forest plots of the results
#' @export 
setGeneric("calc_lmm", function(object, which_data, formula=NULL, cores=NULL) standardGeneric("calc_lmm"))
setMethod("calc_lmm" , "metime_analyser", function(object, which_data, formula=NULL, cores=NULL){
  #define number of cores
  if(is.null(cores)){
    ncores = parallel::detectCores()-1
  }else{
    ncores = as.numeric(cores)
  }
  
  # add sanity checks 
  if(is.null(formula) || !which_data %in% names(object@list_of_data)){
    my_results = data.frame()
    warning("calc_lmm() could not match formula or data")
  }
  
  # add check if all variables are in formula
  
  my_results=lapply(formula, function(x){
    cat(x)
    my_var = x %>% 
      gsub(pattern=" ", replacement="") %>% 
      strsplit(x, split="[+]|~") %>% 
      unlist()
    my_var[length(my_var)] =  my_var[length(my_var)] %>% gsub(pattern="[(1|)]| [)]", replacement="")
    my_data <- object@list_of_data[[which_data]] %>% 
      dplyr::select(all_of(my_var)) %>% 
      na.omit() %>% 
      as.data.frame() 
    
    if(length(unique(my_data$timepoint))==1){
      model_out=data.frame(
        met=my_var[1],
        trait=my_var[2],
        beta = NA,
        pval = NA,
        stringsAsFactors = F
      )
    }else{
      my_model <- lmerTest::lmer(data=my_data,
                                           formula = as.formula(x)
      )
      my_model <- summary(my_model)
      model_out=data.frame(
        met=my_var[1],
        trait=my_var[2],
        beta = my_model[["coefficients"]][2,1],
        pval = my_model[["coefficients"]][2,5],
        stringsAsFactors = F
      )
    }
  })
  my_results_formated = my_results %>% 
    do.call(what=rbind.data.frame) %>% 
    dplyr::mutate(x = beta, 
                  y = met,
                  color=NA,
                  split=trait)
  out <- get_make_plotter_object(
    data = my_results_formated,
    metadata = NULL,
    calc_type="regression",
    calc_info="lmerTest::lmer",
    plot_type="forest",
    style="ggplot"
  )
  return(out)
})

