#' Function to perform matrixEQTL style regression longitudinally 
#' @param object An S4 object of class metime_analyser
#' @param which_data character to define dataset to be used 
#' @param covariates covariates to be used
#' @param feature_selection results from feature_selection
#' @param regression_phenotype character to define which phenotype
#' @param cols_for_meta colnames to be added for meta information. list of character vectors of length which_data
#' @param name character vector for naming the results
#' @param stratifications List to stratify data into a subset. Usage list(name=value)
#' @return plotter object with a forest plot 
#' @export
setGeneric("calc_lm_matrixeqtl", function(object, which_data, covariates, feature_selection, regression_phenotype, name, stratifications, cols_for_meta) standardGeneric("calc_lm_matrixeqtl"))
setMethod("calc_lm_matrixeqtl", "metime_analyser", function(object, which_data, covariates, feature_selection, name, stratifications, cols_for_meta) {
    stopifnot(all(which_data %in% names(object@list_of_data)))
    stopifnot(timepoints %in% names(feature_selection))
    for(i in timepoints) {   
      this_strat_data <- object
      my_adjusted_medication<- feature_selection[[paste0("tp_",i)]]
      my_regression_phenotypes <- this_strat_data@list_of_col_data$regression_data %>% 
            dplyr::filter(type %in% grep(x=type, pattern="regression", value=T)) %>% 
            dplyr::pull(id)
      this_data <- this_strat_data@list_of_data$regression_data %>%
            dplyr::mutate(id=rownames(.[])) %>%
            dplyr::left_join(this_strat_data@list_of_row_data$regression_data[,c("id","time","subject")]) %>% 
            dplyr::filter(time==paste0("t",i))
      if(!class(this_data[ ,regression_phenotype]) %in% c("numeric", "integer")) {
          this_data[ ,regression_phenotype] <- as.factor(this_data[ ,regression_phenotype]) %>% as.numeric
      }

        results_adjusted[[paste0("tp",i)]]  = lapply(unique(this_strat_data@list_of_col_data[[which_data]]$id), function(x) {
          cat(x)
          this_covariates = c(covariates, my_adjusted_medication$id_y[which(my_adjusted_medication$id_x==x)])
    
          this_calc_data <- this_data %>% 
            dplyr::select(all_of(c(x,my_regression_phenotypes,this_covariates)))
          metabolite_data = this_data %>% dplyr::select(x)
          trait_data = this_data %>%
                    dplyr::select(all_of(my_regression_phenotypes))
          covariate_data = this_data %>%
                    dplyr::select(all_of(this_covariates))

          length = intersect(rownames(metabolite_data), rownames(trait_data)) %>% 
                    intersect(rownames(covariate_data)) %>% length()
          real_length = rownames(metabolite_data) %>% length()
          stopifnot(length==real_length)
          met_data <- MatrixEQTL::SlicedData$new()
          met_data$CreateFromMatrix(t(metabolite_data[use_rows,])) # leave in features X metabolites
          t_data <- MatrixEQTL::SlicedData$new()
          t_data$CreateFromMatrix(t(trait_data[use_rows,]))
          cov_data <- MatrixEQTL::SlicedData$new()
          cov_data$CreateFromMatrix(t(covariate_data[use_rows,]))
  
          results_list <- MatrixEQTL::Matrix_eQTL_main(
                snps = met_data, #met_data
                gene = t_data, #t_data
                cvrt = cov_data,
                pvOutputThreshold = 1, # output all p-values
                verbose = F#,
                #useModel=1 # use model 1= modelLINEAR
                )
          out <- results_list$all$eqtls %>%
            dplyr::rename(trait=gene, met=snps) %>%
            dplyr::mutate(met=x,tp=i)
      }) %>% do.call(what=rbind.data.frame)
    }
    all_results <- results_adjusted %>% 
        do.call(what=rbind.data.frame) %>% 
        dplyr::mutate(trait = as.character(trait))
    metadata <- get_metadata_for_columns(object=object, which_data=which_data, columns=cols_for_meta, 
                 names=c("name", "pathway"), index_of_names=rep("id", each=length(which_data)))
    out <- lapply(unique(all_results$trait), function(x) {
            data <- all_results%>% 
               dplyr::filter(pvalue<=0.05/281, 
                              trait==i)
            plotter_object <- get_results(object=object, data=list(data), metadata=list(metadata), 
                                calc_type="linear_model", 
                                calc_info=paste(which_data, "linear_model", "with", 
                                    ifelse(length(stratifications)>=1, 
                                        paste(stratifications, collapse="_"), "full_data"),
                                    sep="_"), 
                    name=name)
      })
    return(out)
    
})