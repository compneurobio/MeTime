#' Comparison of regression analysis
#' @description A function to compare results from regression analysis (LM, LMM, GAMM)
#' @param object a S4 object of class metime_analyser.
#' @param method a character vector of methods 'sign' (sign of beta), 'het' (heterogeneity), 'cor' (correlation). See details for more information. 
#' @param result_index a character vector of results to compare. Default set to NULL. If set to NULL all results from regression are used.
#' @param name a character containing the name of the new result element appended to the analyzer object. Default set to 'comp_regression_1'.
#' @param update_table a logical whether the results should be added to the plotting section. If true tables will be printed in the report. 
#' @details Multiple options are available to compare betas and betas + se from regression analysis. \n sign: comparison of beta signs. \n het (heterogeneity): Significance of heterogeneity across two models as described by \href{<https://doi.org/10.1111/j.1745-9125.1998.tb01268.x>}{Paternoster, 1998} and \href{<https://doi.org/10.1371/journal.pgen.1002215>}{Mittelstrass, 2011} . \n cor (correlation): Similarity of betas based on Pearson correlation
#' @return a S4 object of class metime_analyser object with function output appended to results
#' @seealso [meta_parametric_test], [meta_nonparametric_test], [meta_network]
#' @export  
setGeneric("meta_regression", function(object,  method=c("sign","cor","het"), result_index=NULL, name= "meta_regression_1", update_table = T) standardGeneric("meta_regression"))
setMethod("meta_regression", "metime_analyser", function(object, method=c("sign","cor","het"), result_index=NULL, name= "meta_regression_1", update_table = T){
  if(!all(method %in% c("sign", "cor","het"))) {
    warning('method has to be one of "sign" (sign of beta), "cor" (correlation), and "het" (heterogeneity). Exiting without making any changes')
    return(object)
  }else{
    # save file to be modified
    out <- list()
    
    if(name %in% names(object@results)) {
      warning("name of the results was previously used, using a different name")
      index <- name %>% gsub(pattern="meta_regression_", replacement="") %>% as.numeric() +1
      name <- paste0("meta_regression_",index)
    }
    
    # Preprocessing ----
    ## get all result indices if result_index is null
    if(is.null(result_index)){
      result_index <- lapply(names(object@results), function(x) if(all(object@results[[x]]$information$calc_type =="regression")) x) %>% unlist()
    }
    
    my_data <- sapply(object@results[result_index], "[", "plot_data", simplify = T) %>% unlist(recursive = F)  # compile data for comparison
    names(my_data) <- names(my_data) %>% gsub(pattern=".plot_data", replacement = "") # change model names
    my_combn <- combn(names(my_data), 2) %>% t() %>% as.data.frame() %>% setNames(c("result1","result2")) # get combinations
    
    # Compare by sign ----
    if("sign" %in% method){
      ## calculate result
      this_out <- lapply(1:nrow(my_combn), function(x){
        this_result<-full_join(my_data[[my_combn$result1[x]]] %>%  # join data by met column
                                 dplyr::select(met, trait, beta) %>% 
                                 dplyr::rename(met=met, trait1=trait, beta1=beta) %>% 
                                 dplyr::mutate(sign1=ifelse(beta1>=0, "+","-")),
                               my_data[[my_combn$result2[x]]] %>% 
                                 dplyr::select(met, trait, beta) %>% 
                                 dplyr::rename(met=met, trait2=trait, beta2=beta)%>% 
                                 dplyr::mutate(sign2=ifelse(beta2>=0, "+","-")),
                               by="met") %>% 
          dplyr::mutate(model1=my_combn$result1[x], model2=my_combn$result2[x], combined = paste0(sign1," ", sign2)) %>% # add colums needed for later comparison
          dplyr::count(model1, model2, combined) %>% 
          tidyr::spread(key = combined, value = n)
      }) %>% 
        plyr::rbind.fill()
      out[["sign"]] <- this_out
    }
    
    # Compare by heterogeneity ----
    if("het" %in% method){
      this_out <- lapply(1:nrow(my_combn), function(x){
        if(all(c("se","pval","beta") %in% names(my_data[[my_combn$result1[x]]])) & all(c("se","pval","beta") %in% names(my_data[[my_combn$result2[x]]]))  ){
          #pval and se have to be in the data frames
          full_join(my_data[[my_combn$result1[x]]] %>% 
                      dplyr::select(met, trait, beta, se) %>% 
                      dplyr::rename(met=met, trait1=trait, beta1=beta, se1=se),
                    my_data[[my_combn$result2[x]]] %>% 
                      dplyr::select(met, trait, beta, se) %>% 
                      dplyr::rename(met=met, trait2=trait, beta2=beta, se2=se),
                    by="met") %>% 
            dplyr::mutate(diff.t =  (beta1 - beta2)/sqrt(se1^2 + se2^2)) %>% 
            dplyr::mutate(diff.p  = 2*pnorm(- abs(diff.t)), 
                          i_sq = ifelse(abs(diff.t)>1, ((abs(diff.t) - 1) / abs(diff.t))*100, 0)) %>% 
            dplyr::mutate(significant = ifelse(diff.p<=0.05/nrow(.[]), "non_significant", "significant"))%>%  
            dplyr::count(significant) %>%
            dplyr::mutate(model1=my_combn$result1[x], model2 = my_combn$result2[x])%>% 
            tidyr::spread(key = significant, value = n)
        }else{
        }
      })  %>% 
        plyr::rbind.fill()
      out[["het"]] <- this_out
    }
    
    # compare by correlation ----
    if("cor" %in% method){
      this_out <- lapply(1:nrow(my_combn), function(x){
        
        this_data <- full_join(my_data[[my_combn$result1[x]]] %>% 
                                 dplyr::select(met, trait, beta, se) %>% 
                                 dplyr::rename(met=met, trait1=trait, beta1=beta, se1=se),
                               my_data[[my_combn$result2[x]]] %>% 
                                 dplyr::select(met, trait, beta, se) %>% 
                                 dplyr::rename(met=met, trait2=trait, beta2=beta, se2=se),
                               by="met") 
        
        data.frame(
          model1 = my_combn$result1[x],
          model2 = my_combn$result2[x],
          cor = cor(x=this_data$beta1,y=this_data$beta2, use = "everything",method="pearson"),
          stringsAsFactors = F)
      }) %>% 
        plyr::rbind.fill()
      
      out[["cor"]] <- this_out
    }
    
    # add plot_data ----
    out_object <- get_make_results(object = object, 
                                   data = out, 
                                   metadata = NULL, 
                                   calc_type = rep("meta_regression",length(out)),
                                   calc_info = paste0("meta_regression analysis of",names(out)),
                                   name = name) %>%
      add_function_info(function_name = "meta_regression", 
                        params = list(result_index=result_index,
                                      update_table = update_table,
                                      method=method)
      )
    
    if(update_table){
      out_object@results[[name]][["plots"]][[1]] <- out
    }
  }
    return(out_object)
  }
)
        