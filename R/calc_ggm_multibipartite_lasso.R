
#' An automated fucntion to calculate GGM from multibipartite lasso approach
#' @description automated funtion that can be applied on s4 object of class metime_analyser to calculate a network using
#' multibipartite lasso
#' @param object S4 object of cĺass metime_analyser
#' @param which_data a character or a character vector naming the datasets of interest
#' @param alpha tuning parameter for lasso + ridge regression in glmnet. Default set to 1 to perform LASSO
#' @param nfolds nfolds for cv.glmnet. Default set to 3
#' @param stratifications List to stratify data into a subset. Usage list(name=value)
#' @param cols_for_meta a list of character vectors of column names to be used for visualization of the networks.
#' @param name character to define the name of the results
#' @param cores number of cores to be used in parallel::mclapply(). Default set to 4. 
#' @param ... additional arguments for cv.glmnet function
#' @returns Analyser object with updated results of this calculation 
#' @export
setGeneric("calc_ggm_multibipartite_lasso", function(object, which_data, alpha=1, nfolds=3, stratifications, cols_for_meta, cores=4) standardGeneric("calc_ggm_multibipartite_lasso"))
setMethod("calc_ggm_multibipartite_lasso", "metime_analyser", function(object, which_data, alpha=1, nfolds=3, stratifications, cols_for_meta, cores=4) {
        if(length(which_data)==1) warning("Only one dataset(platform) is being used")

        data_lists <- lapply(which_data, function(a) {
                data_list <- get_stratified_data(object=object, which_data=a,
                    stratifications=stratifications)
                return(data_list)
            })

        all_samples <- lapply(seq_along(data_lists), function(b) {
                samples <- data_lists[[b]][["data"]] %>% rownames() 
                return(samples)
            })

        common_samples <- Reduce(intersect, all_samples)

        data_lists <- lapply(seq_along(data_lists), function(a) {
                data <- data_lists[[a]][["data"]]
                row_data <- data_lists[[a]][["row_data"]]
                data <- data[rownames(data) %in% common_samples, ]
                row_data <- row_data[rownames(row_data) %in% common_samples, ]
                return(list(data=data, row_data=row_data))
            })

        list_of_mats <- lapply(seq_along(data_lists), function(a) {
                  return(data_lists[[a]][["data"]] %>% as.matrix())
            })

        get_betas_for_multibipartite_lasso <- function(list_of_mats, # list of matrices that are divided based on platform or timepoints
           alpha, # alpha parameter for glmnet
           nfolds # nfolds parameter for glmnet
           ) {
          #creating a list to store the data from glmnet
          #code exactly similar to the usual MLP 
          out <- list() # list to store the regression information for each metabolte
          out <- lapply(seq_along(list_of_mats), function(x) {
            each_result <- lapply(seq_along(list_of_mats), function(y) {
                if(x!=y) {
                  x_mat <- as.matrix(list_of_mats[[x]])
                  y_mat <- as.matrix(list_of_mats[[y]])
                  fit_lists <- parallel::mclapply(seq_along(y_mat), function(z) {
                        fit_list <- glmnet::cv.glmnet(x=x_mat, y=y_mat[,z], alpha=alpha, nfolds=nfolds, ...)
                        return(fit_list)
                    }, max.cores=cores)
                  names(fit_lists) <- colnames(y_mat)
                } else {
                  mat <- as.matrix(list_of_mats[[x]])
                  fit_lists <- parallel::mclapply(seq_along(mat), function(z) {
                        y <- as.matrix(mat[,z])
                        x <- as.matrix(mat[,-z])
                        fit_list <- glmnet::cv.glmnet(x=x, y=y, alpha=alpha, nfolds=nfolds, ...)
                        return(fit_list)
                    }, max.cores=cores)
                  names(fit_lists) <- colnames(mat)
                }
                return(fit_lists)
              })
              names(each_result) <- paste(names(list_of_mats)[x], names(list_of_mats), sep="-")
              return(each_result)
            })
          return(out)
        }

        results_list <- get_betas_for_multibipartite_lasso(list_of_mats=list_of_mats, 
          alpha=alpha, nfolds=nfolds, ...)
        
        edge_lists <- lapply(seq_along(results_list), function(a) {
                edge_list <- lapply(seq_along(results_list[[a]]), function(b) {
                      coeffs <- coef(results_list[[a]][[b]])[,1]
                      #coeffs is the vector with all the coefficient values
                      #removing zero coefficients
                      coeffs <- coeffs[!(coeffs==0)]
                      #removing the intercept data
                      coeffs <- coeffs[-1]
                      #source is the metab names
                      node1 <- rep(names(results_list[[a]])[b], each=length(coeffs))
                      #target is the name of the metabolites that have non-zero coeffs
                      node2 <- names(coeffs)
                      # storing it as a list and returning
                      results <- cbind.data.frame(node1, node2, coeffs)
                      return(results)
                  }) %>% do.call(what=rbind.data.frame)
                return(edge_list)
          }) %>% do.call(what=rbind.data.frame)
        edge_lists$node1 <- edge_lists$node1 %>% as.character()
        edge_lists$node2 <- edge_lists$node2 %>% as.character()
        edge_lists$coeffs <- edge_lists$coeffs %>% as.character() %>% as.numeric()
        uids <- apply(edge_lists[ ,c("node1", "node2")], 1, function(x) {
                uid <- paste(stringr::str_sort(x), collapse="_")
                return(uid)
          })
        edge_lists <- edge_lists %>% dplyr::mutate(uid=as.character(uids))
        edge_lists <- reshape(transform(edge_lists, time=ave(coeffs, uids, FUN=seq_along)), 
                      idvar="uids", direction="wide")
        edge_lists$node1 <- gsub(edge_lists$uid, pattern="_[a-z|A-Z|.|0-9]+", replacement="") %>% as.character()
        edge_lists$node2 <- gsub(edge_lists$uid, pattern="[a-z|A-Z|.|0-9]+_", replacement="") %>% as.character()
        metadata <- get_metadata_for_columns(object=object, which_data=which_data, 
            columns=cols_for_meta)
        out <- get_make_results(object=object, data=list(edge_lists), metadata=metadata, 
                calc_type="multibipartite_ggm", calc_info=paste("multibipartite_ggm for ", which_data,
                    "with", ifelse(is.null(stratifications)), "full data", stratifications, sep=" "), 
                name=name, plot_type=list()) %>%
                add_function_info(function_name="calc_ggm_multibipartite_lasso", 
                        params=list(which_data=which_data, stratifications=stratifications, alpha=alpha, nfolds=nfolds, 
                        cols_for_meta=cols_for_meta, name=name))
        return(out)
  })








