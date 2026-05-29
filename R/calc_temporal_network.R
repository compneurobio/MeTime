
#' An automated function to caluclate temporal network with lagged model
#' @description calculates temporal networks for each dataset with a lagged model as used in graphical VAR
#' @param object S4 object of class metime_analyser
#' @param lag which lagged model to use. 1 means one-lagged model, similarly 2,3,..etc
#' @param which_data dataset or datasets to be used
#' @param stratifications List to stratify data into a subset. Usage list(name=value)
#' @param alpha parameter for regression coefficient. Set to 1 for lasso regression
#' @param nfolds nfolds parameter for glmnet style of regression. Default is set to 3
#' @param cols_for_meta a list of character vectors of column names to be used for visualization of the networks.
#' @param names character vector with the same length as that of possible models 
#' @param cluster_profile Character string controlling worker initialization profile.
#'   One of `"auto"`, `"local"`, or `"hpc"`.
#'   - `"auto"`: detects common scheduler environments (e.g. SLURM/PBS/LSF) and chooses `"hpc"` when detected.
#'   - `"local"`: standard local machine setup with default library paths.
#'   - `"hpc"`: cluster-oriented setup; can prepend worker library paths via `hpc_libpaths`.
#'
#' @param hpc_libpaths Optional character vector of library paths to prepend on worker nodes
#'   when `cluster_profile = "hpc"`. Use this when compute nodes need explicit `.libPaths()`
#'   (e.g. non-shared user libraries). Ignored for `"local"` profile.
#' @param num_cores numeric input to define the number of cores that you want to use for parallel computing. Default is set to NULL which is parallel::detectCores() -1.
#' @param ... additional arguments for cv.glmnet function
#' @returns S4 object with updated temporal network results
#' @export
setGeneric("calc_temporal_network", function(object, which_data, lag, stratifications, alpha=1, nfolds=3, cols_for_meta, names, cluster_profile = c("auto", "local", "hpc"), num_cores=NULL,
hpc_libpaths = NULL, ...) standardGeneric("calc_temporal_network"))
setMethod("calc_temporal_network", "metime_analyser", function(object, which_data, lag, stratifications, alpha=1, nfolds=3, cols_for_meta, names, cluster_profile = c("auto", "local", "hpc"), num_cores=NULL,
hpc_libpaths = NULL, ...) {
        cluster_profile = match.arg(cluster_profile)
        if(is.null(stratifications)) {
          times <- object@list_of_row_data[[which_data[1]]]$time %>% unique()
        } else {
          times <- stratifications$time
        }
        data_list <- get_stratified_data(object=object, which_data=which_data,
                    stratifications=stratifications)
        data <- data_list[["data"]]
        row_data <- data_list[["row_data"]]

        data$time <- rownames(data) %>% gsub(pattern="[a-z|A-Z][0-9]+_", replacement="")

        final_data <- lapply(times, function(x) {
                new_data <- data[data$time %in% x, ]
                new_data$time <- NULL
                return(new_data)
            })

        names(final_data) <- times
        count <- 1
        model_seqs <- list()
        while(lag + count <= length(times)) {
            model_seqs[[count]] <- times[seq(from=lag+count, to=count, by=-1)]
            count <- count + 1
        }
        
        names_data <- lapply(model_seqs, function(a) {
                    a <- paste(paste(a, collapse="-"), which_data, sep="_")
                    return(a)
                })

        if(!length(names)==length(model_seqs)) {
            warning("Names are not of correct length and hence using new names")
            names <- names_data
        }
            
        
        models <- list()
        out <- lapply(seq_along(model_seqs), function(x) {
                   network_data <- final_data[names(final_data) %in% model_seqs[[x]]]
                   count <- 1
                   network_data <- lapply(network_data, function(y) {
                           colnames(y) <- paste(colnames(y), model_seqs[[x]][count], sep="_time:")
                           rownames(y) <- rownames(y) %>% gsub(pattern="_[a-z|A-Z][0-9]+", replacement="")
                           count <<- count +1 
                           return(y) 
                     })
                   list_of_names <- lapply(network_data, function(y) {
                       return(rownames(y))
                   })
                   common_samples <- Reduce(intersect, list_of_names)
                   network_data <- lapply(network_data, function(y) {
                       y <- y[rownames(y) %in% common_samples, ]
                       y <- y[order(rownames(y)), ]
                       return(y)
                   })
                   network_data <- do.call(cbind, unname(network_data))
                   ymat <- network_data[ ,grep(model_seqs[[x]][1] ,colnames(network_data), value=TRUE)]
                   xmat <- network_data[ ,!(colnames(network_data) %in% colnames(ymat))]
                   ymat <- as.matrix(na.omit(ymat))
                   xmat <- as.matrix(na.omit(xmat))
                   if(is.null(num_cores)) num_cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)-1
                   cl <- .init_cluster(
                        num_cores = num_cores,
                        profile = cluster_profile,
                        hpc_libpaths = hpc_libpaths,
                        worker_packages = c("glmnet"),
                        export_vars = c("xmat", "ymat", "alpha", "nfolds"),
                        export_env = environment()
                    )
                    on.exit(.stop_cluster(cl), add = TRUE)

                    targets <- colnames(ymat)
                    results <- .apply_with_progress(targets, cl = cl, FUN = function(y) {
                            fit <- glmnet::cv.glmnet(x = xmat, y = ymat[, y], alpha = alpha, nfolds = nfolds, ...)
                            coeffs <- coef(fit)[, 1]
                            coeffs <- coeffs[coeffs != 0][-1]
                            if (length(coeffs) == 0) return(NULL)
                            data.frame(from = names(coeffs), to = rep(y, length(coeffs)), coeffs = as.numeric(coeffs))
                        })
                    results <- do.call(rbind.data.frame, Filter(Negate(is.null), results))
                   results$label <- paste(unlist(lapply(strsplit(as.character(results$from), split="_time:"), function(x) return(x[2]))),
                                        unlist(lapply(strsplit(as.character(results$to), split="_time:"), function(x) return(x[2]))),
                                        results$coeffs, 
                                        sep="-")
                   results$from <- unlist(lapply(strsplit(as.character(results$from), split="_time:"), function(x) return(x[1])))
                   results$to <- unlist(lapply(strsplit(as.character(results$to), split="_time:"), function(x) return(x[1])))  
                   return(results)
            })
        metadata <- get_metadata_for_columns(object=object, which_data=which_data, 
            columns=cols_for_meta)
        names(out) <- names(data)
        for(i in seq_along(model_seqs)) {
            data <- list(out[[i]])
            names(data) <- names_data[i]
            object <- get_make_results(object=object, data=data, metadata=metadata, 
                calc_type="temporal_network", calc_info=paste("temporal_network", model_seqs[[i]], "for", which_data,
                    "with", ifelse(is.null(stratifications)), "full data", stratifications, sep=" "), 
                name=names[i])
            object <- add_function_info(object=object, function_name="calc_temporal_network", 
                        params=list(which_data=which_data, lag=lag, 
                        stratifications=stratifications, alpha=alpha, nfolds=nfolds, 
                        cols_for_meta=cols_for_meta, cores=cores, names=names[i], ...))
        }
        out <- object
        return(out)
  })

