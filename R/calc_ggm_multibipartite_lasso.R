
#' An automated fucntion to calculate GGM from multibipartite lasso approach
#' @description automated funtion that can be applied on s4 object of class metime_analyser to calculate a network using
#' multibipartite lasso
#' @param object S4 object of cĺass metime_analyser
#' @param which_data a character or a character vector naming the datasets of interest
#' @param alpha tuning parameter for lasso + ridge regression in glmnet
#' @param nfolds nfolds for cv.glmnet
#' @param timepoints timepoints of interest that are to be used to build networks(as per timepoints in rows)
#' @param cols_for_meta a list of character vectors of column names to be used for visualization of the networks.
#' @return list of plotter objects that can be used for plotting.
#' @export
setGeneric("calc_ggm_multibipartite_lasso", function(object, which_data, alpha, nfolds, timepoints, cols_for_meta) standardGeneric("calc_ggm_multibipartite_lasso"))
setMethod("calc_ggm_multibipartite_lasso", "metime_analyser", function(object, which_data, alpha, nfolds, timepoints, cols_for_meta) {
          if(length(which_data) > 1) object <- mod_extract_common_samples(object)
          object <- mod_split_acc_to_time(object)
          converted <- mod_code_metab_names(object=object, which_data=which_data)
          new_object <- converted$object
          table <- converted$table
          list_of_data <- new_object@list_of_data[names(list_of_data) %in% which_data]
          list_of_data <- lapply(list_of_data, function(x) return(x[names(x) %in% timepoints]))
          count <- 1
          edge_lists <- list()
          for(i in 1:length(timepoints)) {
              list_of_mats <- list()
              list_of_mats <- lapply(list_of_data, function(x) {
                            x <- x[names(x) %in% timepoints[i]]
                            x <- lapply(x, function(a) return(a[order(rownames(a)), ]))
                            x <- lapply(x, function(a) return(as.matrix(a)))
                            return(x[[1]]) 
              })
              names(list_of_mats) <- which_data
              #glmnet results are stored here in a list
              list <- get_betas_for_multibipartite_lasso(list_of_mats=list_of_mats, 
                                            alpha=alpha, nfolds=nfolds)
              #list to check and collect all the edges that have non-zero coefficient 
              list_to_check <- list()
              #list here is the list with all the details from glmnet
              for(j in 1:length(list)) {
                  #creating a list to store the edges
                  list_of_edges <- list()
                  #list[[1]] is the data for first metabolite
                  for(k in 1:length(list[[j]])) {
                      #coeffs is the vector with all the coefficient values
                      coeffs <- coef(list[[j]][[k]])[,1]
                      #removing zero coefficients
                      coeffs <- coeffs[!(coeffs==0)]
                      #removing the intercept data
                      coeffs <- coeffs[-1]
                      #source is the metab names
                      source <- rep(names(list[[j]])[k], each=length(coeffs))
                      #target is the name of the metabolites that have non-zero coeffs
                      target <- names(coeffs)
                      #storing data as a nested list
                      list_of_edges[[k]] <- cbind(source, target, coeffs) 
                      names(list_of_edges)[k] <- source[1]
                     
                }
                #list_to_check will now store the edges
                list_to_check[[j]] <- list_of_edges
              }
            edge_list <- unlist(list_to_check, recursive=FALSE)
            edge_list <- as.data.frame(do.call(rbind, edge_list)) 
            uids <- c()
            #ERROR FOUND HERE AND CORRECTED!!!!!!!!
            for(m in 1:length(edge_list[,1])) {
                vec <- c()
                vec[1] <- as.character(edge_list$source[m])
                vec[2] <- as.character(edge_list$target[m])
                vec <- vec[order(vec)]
                uids[m] <- paste(vec, collapse="_")
            }
            edge_list <- cbind(edge_list, uids)
            edge_list <- as.data.frame(edge_list)
            colnames(edge_list) <- c("node1", "node2", "coeffs", "uids")
            check <- edge_list[,c("uids", "coeffs")]
            check$uids <- as.character(check$uids)
            check$coeffs <- as.numeric(as.character(check$coeffs))
            check <- reshape(transform(check, time=ave(coeffs, uids, FUN=seq_along)), idvar="uids", direction="wide")
            final_edge_list <- na.omit(check)
            final_edge_list$uids <- as.character(final_edge_list$uids)
            final_edge_list$node1 <- unlist(lapply(strsplit(final_edge_list$uids, split="_"), function(x) return(x[1])))
            final_edge_list$node2 <- unlist(lapply(strsplit(final_edge_list$uids, split="_"), function(x) return(x[2])))
            final_edge_list$uids <- NULL
            edge_lists[[i]] <- final_edge_list            
          }
          names(edge_lists) <- timepoints
          metadata <- get_metadata_for_columns(object=object, which_data=which_data, columns=cols_for_meta, 
                 names=c("name", "pathway"), index_of_names=rep("id", each=length(which_data)))
          out <- lapply(edge_lists, function(x) {
                for(i in 1:length(x$node1)) {
                  x$node1[i] <- table$metabolite[table$id %in% x$node1[i]]
                  x$node2[i] <- table$metabolite[table$id %in% x$node2[i]]
                }
                plotter_object <- get_make_plotter_object(data=x, metadata=metadata, calc_type="multibipartite_ggm", 
                  calc_info=paste("multibipartite lasso network for:", paste(which_data, collapse=" & "), sep=" "), 
                plot_type="network", style="visNetwork")
                return(plotter_object) 
            })
          return(out)
  })


