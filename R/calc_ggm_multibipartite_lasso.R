
#' An automated fucntion to calculate GGM from multibipartite lasso approach
#' @description automated funtion that can be applied on s4 object of class metime_analyser to calculate a network using
#' multibipartite lasso
#' Conceptualization & Methodology: Richa Batra and Jan Krumsiek.
#' Implementation: Richa Batra, Bharadwaj Marella, Patrick Weinisch and Matthias Arnold.
#' @param object S4 object of cĺass metime_analyser
#' @param which_data a character or a character vector naming the datasets of interest
#' @param alpha tuning parameter for lasso + ridge regression in glmnet. Default set to 1 to perform LASSO
#' @param nfolds nfolds for cv.glmnet. Default set to 3
#' @param stratifications List to stratify data into a subset. Usage list(name=value)
#' @param cols_for_meta a list of character vectors of column names to be used for visualization of the networks.
#' @param name character to define the name of the results
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
#' @returns Analyser object with updated results of this calculation 
#' @export
setGeneric("calc_ggm_multibipartite_lasso", function(object, which_data, alpha=1, nfolds=3, stratifications, cols_for_meta, cluster_profile = c("auto", "local", "hpc"), num_cores=NULL,
hpc_libpaths = NULL, ...) standardGeneric("calc_ggm_multibipartite_lasso"))
setMethod("calc_ggm_multibipartite_lasso", "metime_analyser", function(object, which_data, alpha=1, nfolds=3, stratifications, cols_for_meta, cluster_profile = c("auto", "local", "hpc"), num_cores=NULL,
hpc_libpaths = NULL, ...) {
        if(length(which_data)==1) warning("Only one dataset(platform) is being used")
        cluster_profile=match.arg(cluster_profile)
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

      get_betas_for_multibipartite_lasso <- function(
          list_of_mats,
          alpha,
          nfolds,
          num_cores = NULL,
          cluster_profile = cluster_profile,
          hpc_libpaths = hpc_libpaths,
          ...
        ) {
          cluster_profile <- match.arg(cluster_profile)

          # --- 1) Build global task table: one glmnet fit per row ---
          pair_grid <- expand.grid(
            x = seq_along(list_of_mats),
            y = seq_along(list_of_mats),
            stringsAsFactors = FALSE
          )

          task_list <- lapply(seq_len(nrow(pair_grid)), function(i) {
            x <- pair_grid$x[i]
            y <- pair_grid$y[i]

            # target dimension depends on branch
            n_targets <- if (x != y) {
              ncol(list_of_mats[[y]])
            } else {
              ncol(list_of_mats[[x]])
            }

            data.frame(
              pair_id = i,
              x = x,
              y = y,
              z = seq_len(n_targets),
              stringsAsFactors = FALSE
            )
          })

          tasks <- do.call(rbind, task_list)

          # --- 2) Single cluster for all tasks ---
          cl <- .init_cluster(
            num_cores = num_cores,
            profile = cluster_profile,
            hpc_libpaths = hpc_libpaths,
            worker_packages = c("glmnet"),
            export_vars = c("list_of_mats", "tasks", "alpha", "nfolds"),
            export_env = environment()
          )
          on.exit(.stop_cluster(cl), add = TRUE)

          # --- 3) Parallel over flattened tasks ---
          fit_out <- .apply_with_progress(
            X = seq_len(nrow(tasks)),
            cl = cl,
            FUN = function(i) {
              tt <- tasks[i, ]

              x <- tt$x
              y <- tt$y
              z <- tt$z

              if (x != y) {
                x_mat <- as.matrix(list_of_mats[[x]])
                y_mat <- as.matrix(list_of_mats[[y]])
                y_vec <- y_mat[, z]
                target_name <- colnames(y_mat)[z]
              } else {
                mat <- as.matrix(list_of_mats[[x]])
                y_vec <- mat[, z]
                x_mat <- mat[, -z, drop = FALSE]
                target_name <- colnames(mat)[z]
              }

              fit <- glmnet::cv.glmnet(
                x = x_mat,
                y = y_vec,
                alpha = alpha,
                nfolds = nfolds,
                ...
              )

              list(x = x, y = y, z = z, target_name = target_name, fit = fit)
            }
          )

          # --- 4) Rebuild original nested structure: out[[x]][[y]] -> named list of fits ---
          n <- length(list_of_mats)
          out <- vector("list", n)

          for (x in seq_len(n)) {
            each_result <- vector("list", n)

            for (y in seq_len(n)) {
              sel <- vapply(fit_out, function(a) a$x == x && a$y == y, logical(1))
              cell <- fit_out[sel]

              # keep deterministic order by z
              ord <- order(vapply(cell, `[[`, numeric(1), "z"))
              cell <- cell[ord]

              fits <- lapply(cell, `[[`, "fit")
              names(fits) <- vapply(cell, `[[`, character(1), "target_name")

              each_result[[y]] <- fits
            }

            names(each_result) <- paste(names(list_of_mats)[x], names(list_of_mats), sep = "-")
            out[[x]] <- each_result
          }

          out
        }

        results_list <- get_betas_for_multibipartite_lasso(list_of_mats=list_of_mats, 
          alpha=alpha, nfolds=nfolds,num_cores = num_cores,
          cluster_profile = cluster_profile,
          hpc_libpaths = hpc_libpaths, ...)
        
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
        out <- get_make_results(object=object, data=list(mbpl=edge_lists), metadata=metadata, 
                calc_type="multibipartite_ggm", calc_info=paste("multibipartite_ggm for ", which_data,
                    "with", ifelse(is.null(stratifications)), "full data", stratifications, sep=" "), 
                name=name, plot_type=list()) %>%
                add_function_info(function_name="calc_ggm_multibipartite_lasso", 
                        params=list(which_data=which_data, stratifications=stratifications, alpha=alpha, nfolds=nfolds, 
                        cols_for_meta=cols_for_meta, name=name))
        return(out)
  })








