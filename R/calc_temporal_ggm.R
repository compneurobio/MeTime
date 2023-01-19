

#' An automated function to caluclate temporal network with lagged model
#' @description calculates temporal networks for each dataset with a lagged model as used in graphical VAR
#' @param object S4 object of class metab_analyser
#' @param lag which lagged model to use. 1 means one-lagged model, similary 2,3,..etc
#' @param which_data dataset or datasets to be used
#' @param timepoints timepoints of interest that are to be used to build networks(in the order of measurement)
#' @param alpha parameter for regression coefficient
#' @param nfolds nfolds parameter for glmnet style of regression
#' @param cols_for_meta a list of character vectors of column names to be used for visualization of the networks.
#' @param cores Number of cores to be used for the process
#' @return temporal network data with edgelist and regression values
#' @export
setGeneric("calc_temporal_ggm", function(object, which_data, lag, timepoints, alpha, nfolds, cols_for_meta, cores) standardGeneric("calc_temporal_ggm"))
setMethod("calc_temporal_ggm", "metime_analyser", function(object, which_data, lag, timepoints, alpha, nfolds, cols_for_meta, cores) {
        if(length(which_data) > 1) object <- mod_extract_common_samples(object)
        object <- mod_split_acc_to_time(object)
        list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
        list_of_data <- lapply(list_of_data, function(x) {
              x <- x[names(x) %in% timepoints]
              x <- lapply(x, function(a) return(a[order(rownames(a)), ]))
              return(x)
          })
        final_data <- list()
        list_of_data <- unname(list_of_data)
        for(i in 1:length(timepoints)) {
            timepoint_list <- lapply(list_of_data, function(x) return(x[[i]]))
            final_data[[i]] <- do.call(cbind, timepoint_list)
        }
        names(final_data) <- timepoints
        count <- 1
        model_seqs <- list()
        while(lag + count <= length(timepoints)) {
            model_seqs[[count]] <- timepoints[seq(from=lag+count, to=count, by=-1)]
            count <- count + 1
        }
        models <- list()
        out <- lapply(seq_along(model_seqs), function(x) {
                   network_data <- final_data[names(final_data) %in% model_seqs[[x]]]
                   count <- 1
                   network_data <- lapply(network_data, function(y) {
                           colnames(y) <- paste(colnames(y), model_seqs[[x]][count], sep="_time:")
                           rownames(y) <- unlist(lapply(strsplit(rownames(x), split="_"), function(x) return(x[1])))
                           count <<- count +1 
                           return(x) 
                     })
                   network_data <- unname(network_data)
                   list_of_names <- lapply(network_data, function(y) {
                       return(rownames(y))
                   })
                   common_samples <- Reduce(intersect, list_of_names)
                   network_data <- lapply(network_data, function(y) {
                       y <- y[rownames(y) %in% common_samples, ]
                       y <- y[order(rownames(y)), ]
                       return(y)
                   })
                   network_data <- do.call(cbind, network_data)
                   ymat <- network_data[ ,grep(model_seqs[[x]][1] ,colnames(network_data), value=TRUE)]
                   xmat <- network_data[ ,!(colnames(network_data) %in% colnames(ymat))]
                   ymat <- as.matrix(na.omit(ymat))
                   xmat <- as.matrix(na.omit(xmat))
                   results <- parallel::mclapply(seq_along(colnames(a)), function(y) {
                            result <- glmnet::cv.glmnet(x=xmat, y=ymat[,y], alpha=alpha, nfolds=nfolds)
                            coeffs <- coef(result)[,1]
                            coeffs <- coeffs[!(coeffs==0)]
                            coeffs <- coeffs[-1]
                            target <- names(coeffs)
                            source <- rep(colnames(ymat)[y], each=length(coeffs))
                            result <- as.data.frame(cbind(source, target, coeffs))
                            return(result)
                    }, max.cores=cores) %>% do.call(what=rbind.data.frame)
                   colnames(results) <- c("node1", "node2", "coeffs")
                   results$label <- paste(unlist(lapply(strsplit(as.character(results$node1), split="_time:"), function(x) return(x[2]))),
                                              unlist(lapply(strsplit(as.character(results$node2), split="_time:"), function(x) return(x[2]))),
                                                sep="-")
                   results$node1 <- unlist(lapply(strsplit(as.character(results$node1), split="_time:"), function(x) return(x[1])))
                   results[[i]]$node2 <- unlist(lapply(strsplit(as.character(results$node2), split="_time:"), function(x) return(x[1])))  
                   return(results)
            })
        models <- lapply(seq_along(model_seqs), function(x) {
                    names(out)[x] <- paste(model_seqs[[x]], collapse="-")
                    return(out)
            })
        metadata <- get_metadata_for_columns(object=object, which_data=which_data, columns=cols_for_meta, 
                 names=c("name", "pathway"), index_of_names=rep("id", each=length(which_data)))
        out <- lapply(models, function(x) {
              plotter_object <- get_make_plotter_object(data=x, metadata=metadata, calc_type="temporal_network", 
                  calc_info=paste("Temporal network for:", paste(which_data, collapse=" & "), sep=" "), 
                plot_type="network", style="visNetwork", aesthetics=NULL)
              return(plotter_object)
          })
        return(out)
  })


