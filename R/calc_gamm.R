
#' Function to perform Generalized additive models
#' @description Function to perform Generalized additive models
#' @param object An S4 object of class metime_analyser
#' @param which_data Dataset to be used for this analysis
#' @param formula A dataframe listing the formulae of the gamms to be used
#' @param cores number of cores of the system to be used. Can also be set to NULL
#' @return plotter object with GAM results
#' @export
setGeneric("calc_gamm", function(object, which_data, formula, cores) standardGeneric("calc_gamm"))
setMethod("calc_gamm", "metime_analyser", function(object, which_data, formula, cores) {
  if(is.null(cores)) {
    ncores=parallel::detectCores()-1
  } else {
    ncores=as.numeric(cores)
  }
  
  # add sanity checks 
  if(is.null(formula) || !which_data %in% names(object@list_of_data)) {
    my_results = data.frame()
    warning("calc_gamm() could not match formula or data")
  }
  
  # add check if all variables are in formula
  my_results <- parallel::mclapply(formula, mc.cores=ncores, function(x) {
      cat(which(formula==x), "; ")
      my_var = x %>% 
      gsub(pattern=" ", replacement="") %>% 
      gsub(pattern="s(", replacement="", fixed = T) %>% 
      strsplit(x, split="[+]|~") %>% 
      unlist() %>% 
      gsub(pattern=")", replacement="", fixed = T) %>%
      gsub(pattern=',bs=\"re\"', replacement="", fixed = T)
    
      my_data <- object@list_of_data[[which_data]] %>% 
      dplyr::select(all_of(my_var))
    
      if(length(unique(my_data$time)) == 1) {
        model_out <- data.frame(
          met = my_var[1],
          trait = my_var[2],
          beta = NA,
          pval = NA,
          stringsAsFactors = F)
      } else {
        this_model <-try(mgcv::gamm(data = my_data,
                              formula= as.formula(x), 
                              method="REML"),
                          silent = T)
        if(class(this_model) == "try-error") {
          model_out <- data.frame(
            met = my_var[1],
            trait = my_var[2],
            beta = NA,
            pval = NA,
            stringsAsFactors = F)
        } else {
          my_model <- summary(this_model$gam)
          model_out = data.frame(
            met=my_var[1],
            trait=my_var[2],
            beta = my_model[["p.coeff"]][[2]],
            pval = my_model[["p.pv"]][[2]],
            stringsAsFactors = F)
        }
      }
      return(model_out)
  })
  my_results_formated <- my_results %>% 
      do.call(what=rbind.data.frame) %>% 
      dplyr::mutate(x = beta, 
                  y = met,
                  color=NA,
                  split=trait) %>% 
      dplyr::select(x, y, color, split)
  rownames(my_results_formated) = 1:nrow(my_results_formated)
  
  out <- get_make_plotter_object(
    data = my_results_formated,
    metadata = NULL,
    calc_type="regression",
    calc_info="mgcv::gamm",
    plot_type="forest",
    style="ggplot"
  )
  return(out)
})

