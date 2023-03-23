#' Function to add the clusters obtained from wgcna
#' @description A function to add cluster assiginment to the col_data from WGCNA
#' @param object An S4 object of class metime_analyser
#' @param which_data character to define which dataset is to be used
#' @param append logical if set to true adds the new data to the object used else creates new object
#' @param results_index results_index if clusters were previously calculated else set to NULL which is also default
#' @param cols_for_meta A list of named character vector to extract col_data to add it for the eigendata
#' @param ... arguments for add_clusters_wgcna
#' @return metime_analyser object with new dataset with eigendata of the metabolites
#' @export  
setGeneric("mod_trans_eigendata", function(object, which_data, append, results_index, cols_for_meta, ...) standardGeneric("mod_trans_eigendata"))
setMethod("mod_trans_eigendata", "metime_analyser", function(object, which_data, append, results_index, cols_for_meta, ...) {
			if(!which_data %in% names(object@list_of_data)) {
          warning("which_data is not found in the datasets of the object. Exiting without making changes")
          return(object)
      }
      if(is.null(results_index)) {
          object <- calc_clusters_wgcna(object=object, which_data=which_data, ...)
          results_index <- length(object@results)
      }
			cluster_info <- object@results[[results_index]][ ,c("id", "modules")]
      colnames(cluster_info)[2] <- "cluster"
			data_list <- split(cluster_info, f = cluster_info$cluster)
			newdata <- object@list_of_data[[which_data]]
        	data_to_append <- lapply(names(data_list), function(x) {
                  modules <- data_list[[x]]
                  dummy <- newdata[ ,modules[,1]]
                  if(!x %in% "0") {
                      pca <- prcomp(dummy, scale.=T, center=T)
                      dummy <- as.data.frame(pca$x[,1])
                      colnames(dummy) <- c(paste("eigenmetabolite_cluster_", x, sep=""))
                  }
                  return(dummy)
          	}) %>% do.call(what=cbind.data.frame)


        	col_info <- lapply(names(data_list), function(x) {
                  modules <- data_list[[x]]
                  if(x %in% "0") {
                      col_info <- data.frame(id=modules[,1], modules=paste("###", modules[,1], sep=""))
                  } else {
                  	  id <- paste("eigenmetabolite_cluster_", x, sep="")
                  	  modules <- paste("###", modules[,1], sep="")
                  	  col_info <- data.frame(id=id, modules=modules)
                  }
                  return(col_info)
          	}) %>% do.call(what=rbind.data.frame)           
            if(is.null(cols_for_meta)) {
                col_data <- object@list_of_col_data[[which_data]]$id %>% as.data.frame()
            } else {
                col_data <- get_metadata_for_columns(object=object, which_data=which_data, 
                     columns=cols_for_meta)
            }
          	col_data_true <- col_data[col_data$id %in% colnames(data_to_append), ]
          	col_data_true <- col_data_true[order(col_data_true$id), ]

          	cluster_data_true <- colnames(data_to_append)[grep("eigen*", colnames(data_to_append))]
          	cluster_data_true <- data.frame(id=cluster_data_true)
          	rownames(cluster_data_true) <- cluster_data_true$id
          
          	col_data_true <- plyr::rbind.fill(col_data_true, cluster_data_true)
          	col_data_true <- col_data_true[order(col_data_true$id), ]
          	col_info <- col_info[order(col_info$id), ]
          	col_data_true <- cbind.data.frame(col_data_true, col_info$modules)
            colnames(col_data_true)[length(colnames(col_data_true))] <- "modules"
          
          	if(append) {
          		out <- add_dataset(object, data=final_data_true, col_data=col_data_true, 
                    row_data=object@list_of_row_data[[which_data]], name=paste(which_data, "with_eigenmetabs", sep="_"))
          	} else {
          		out <- get_make_analyser_object(data=data_to_append, col_data=col_data_true, 
          				row_data=object@list_of_row_data[[which_data]], 
          				annotations_index=list(), name=paste(which_data, "with_eigenmetabs", sep="_"), results=list())
          	}
            out <- add_function_info(object=out, function_name="mod_trans_eigendata", 
                params=list(which_data=which_data, append=append, clusters=clusters, ...))
          	return(out)
	})