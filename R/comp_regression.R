#' Comparison of regression analysis
#' @description A function to compare results from regression analysis (LM, LMM, GAMM)
#' @param object a S4 object of class metime_analyser.
#' @param result_index a character vector of results to compare.
#' @param method a character vector of methods to be applied for comparison. Options: 'sign', 'heterogeneity'. See details for more information.
#' @details Multiple options are available to compare betas and betas + se from regression analysis. Sign: comparison of beta signs. heterogeneity: comparison of betas and se.
#' @return a S4 object of class metime_analyser object with function output appended to results
#' @export  
setGeneric("comp_regression", function(object, result_index=NULL, method="sign") standardGeneric("comp_regression"))
setMethod("comp_regression", "metime_analyser", function(object, result_index=NULL, method="sign"){
  if(length(method) > 1){
    warning("comp_regression() can only use one method")
    out <- object
  }else{
    #save output file
    out <- object
    
    # Preprocessing ----
    ## get all result indices
    if(is.null(result_index)){
      result_index <- lapply(names(out@results), function(x) if(all(out@results[[x]]$information$calc_type =="regression")) x) %>% unlist()
    }
    
    ## compile data for comparison
    my_data <-list()
    for(i in result_index){
      my_data <- append(my_data,out@results[[i]]$plot_data %>% magrittr::set_names(paste0(i,"_",names(.[]))))
    }
    
    ## get data combinations
    my_combn <- combn(names(my_data), 2) %>% t() %>% as.data.frame() %>% 
      magrittr::set_colnames(c("result1","result2"))
    
    # Compare by sign ----
    if("sign" %in% method){
      # if only testing for sign
      this_out <- lapply(1:nrow(my_combn), function(x){
        this_result<-full_join(my_data[[my_combn$result1[x]]] %>% 
                                 dplyr::select(met, trait, beta) %>% 
                                 dplyr::rename(met=met, trait1=trait, beta1=beta) %>% 
                                 dplyr::mutate(sign1=ifelse(beta1>=0, "+","-")),
                               my_data[[my_combn$result2[x]]] %>% 
                                 dplyr::select(met, trait, beta) %>% 
                                 dplyr::rename(met=met, trait2=trait, beta2=beta)%>% 
                                 dplyr::mutate(sign2=ifelse(beta2>=0, "+","-")),
                               by="met") 
        this_result
      })
      names(this_out) <- paste0(my_combn$result1, "_", my_combn$result2)
      ### get name of compariosn file
      name <- "comp_regression_sign_1"
      
      if(name %in% names(out@results)) {
        warning("name of the results was previously used, using a different name")
        index <- name %>% gsub(pattern="comp_regression_sign_", replacement="") %>% as.numeric() +1
        name <- paste0("comp_regression_sign_",index)
      }
      out <- get_make_results(object = out, 
                              data = this_out, 
                              metadata = NULL, 
                              calc_type = rep("table",each=length(this_out)),
                              calc_info = paste0("compare sign of", names(this_out)),
                              name = name) %>%
        add_function_info(function_name = "comp_regression", 
                          params = list(result_index=result_index,
                                        method=method)
        )
      out@results[[name]][["plots"]] <- list(
        list(
          comp_regression_sign=lapply(1:nrow(my_combn), function(x) this_out[[x]] %>% dplyr::count(sign1, sign2) %>% magrittr::set_colnames(c(my_combn[x,],"n")))
        )
      )
      names(out@results[[name]][["plots"]][[1]][[1]]) <- names(this_out)
    }
    # Compare by heterogeneity ----
    if("heterogeneity" %in% method){
      this_out <- lapply(1:nrow(my_combn), function(x){
        if(all(c("se","pval") %in% names(my_data[[my_combn$result1[x]]])) & all(c("se","pval") %in% names(my_data[[my_combn$result2[x]]]))  ){
          #pval and se have to be in the data frames
          this_data<-full_join(my_data[[my_combn$result1[x]]] %>% 
                                 dplyr::select(met, trait, beta, se) %>% 
                                 dplyr::rename(met=met, trait1=trait, beta1=beta, se1=se),
                               my_data[[my_combn$result2[x]]] %>% 
                                 dplyr::select(met, trait, beta, se) %>% 
                                 dplyr::rename(met=met, trait2=trait, beta2=beta, se2=se),
                               by="met")
          
          diff.t <- (this_data$beta1 - this_data$beta2)/sqrt(this_data$se1^2 + this_data$se2^2)
          this_result <- data.frame(
            stringsAsFactors =F,
            met = this_data$met,
            trait_model1 = this_data$trait1,
            trait_model2 = this_data$trait2,
            diff.t = ((this_data$beta1 - this_data$beta2)/sqrt(this_data$se1^2 + this_data$se2^2))) %>% 
            dplyr::mutate(diff.p  = 2*pnorm(- abs(diff.t)), 
                          i_sq = ifelse(abs(diff.t)>1, ((abs(diff.t) - 1) / abs(diff.t))*100, 0)) %>% 
            dplyr::mutate(significant = ifelse(diff.p<=0.05/nrow(.[]), FALSE, TRUE))
          this_result
          # summary <- this_result %>% 
          #   count(significant) %>% 
          #   dplyr::mutate(model1 = my_combn$result1[[x]],
          #                 model2 = my_combn$result2[[x]])
        }else{
          # summary <- data.frame(
          #   stringsAsFactors = F,
          #   significant=NA,
          #   n=NA,
          #   model1=my_combn$result1[[x]],
          #   model2=my_combn$result2[[x]]
          # )
        }
        this_result
      })
      names(this_out) <- paste0(my_combn$result1, "_", my_combn$result2)
      ### get name of compariosn file
      name <- "comp_regression_heterogeneity_1"
      if(name %in% names(out@results)){
        warning("name of the results was previously used, using a different name")
        index <- name %>% gsub(pattern="comp_regression_heterogeneity_", replacement="") %>% as.numeric()+1
        name <- paste0("comp_regression_heterogeneity_",index)
      }
      
      out <- get_make_results(object = out, 
                              data = this_out, 
                              metadata = NULL, 
                              calc_type = rep("table",each=length(this_out)),
                              calc_info = paste0("compare heterogeneity of", names(this_out)),
                              name = name) %>%
        add_function_info(function_name = "comp_regression", 
                          params = list(result_index=result_index,
                                        method=method)
        )
      # add table to plots
      out@results[[name]][["plots"]] <- list(
        list(
          comp_regression_hetero=lapply(1:nrow(my_combn), function(x) this_out[[x]] %>%  
                                        count(significant) %>%
                                        dplyr::mutate(model1=my_combn$result1[x], model2 = my_combn$result2[x])
          )
        )
        )
      names(out@results[[name]][["plots"]][[1]][[1]]) <- names(this_out)
    }
  }
  return(out)
})
        