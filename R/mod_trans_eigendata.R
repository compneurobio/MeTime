#' Function to add the clusters obtained from wgcna
#' @description Modification(mod) function to get cluster assiginment and generate eigenmetabolite matrix from a dataset.
#' @param object An S4 object of class metime_analyser
#' @param which_data character to define which dataset is to be used
#' @param append logical if set to true adds the new data to the object used else creates new object
#' @param results_index results_index if clusters were previously calculated else set to NULL(default)
#' @param cols_for_meta A list of named character vector to extract col_data to add it for the eigendata. will be parsed to 
#' get_metadata_for_columns
#' @seealso [get_metadata_for_columns], [calc_clusters_wgcna]
#' @param ... arguments for calc_clusters_wgcna. Make sure to set the correct baseline value if you are using the function directly
#' @return metime_analyser object with new dataset with eigendata of the metabolites
#' @export  
setGeneric("mod_trans_eigendata", function(object, which_data, append, results_index=NULL, cols_for_meta=NULL, name="WGCNA_clusters_1", ...) standardGeneric("mod_trans_eigendata"))
setMethod("mod_trans_eigendata", "metime_analyser", function(object, which_data, append, results_index=NULL, cols_for_meta=NULL, name="WGCNA_clusters_1", ...) {
	  if(!length(which_data)==1) {
        warning("mod_trans_eigendata(): This function calculates clusters for only one dataset at a time. Exiting without making any changes.")
        return(object)
      }
      if(!which_data %in% names(object@list_of_data)) {
          warning("mod_trans_eigendata(): which_data is not found in the datasets of the object. Exiting without making changes")
          return(object)
      }
      if(is.null(results_index)) {
          object <- calc_clusters_wgcna(object=object, which_data=which_data, cols_for_meta=cols_for_meta, baseline="t0", name=name, ...)
          results_index <- length(object@results)
      }
	  cluster_info <- object@results[[results_index]]$plot_data[[1]][ ,c("id", "modules")]
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
                  	  modules <- paste(c("###", modules[,1]), collapse="###")
                  	  col_info <- data.frame(id=id, modules=modules)
                  }
                  return(col_info)
          	}) %>% do.call(what=rbind.data.frame)           
            if(is.null(cols_for_meta)) {
                col_data <- data.frame(id=object@list_of_col_data[[which_data]]$id)
            } else {
                col_data <- get_metadata_for_columns(object=object, which_data=which_data, 
                     columns=cols_for_meta)
            }
          	col_data_true <- col_data[col_data$id %in% colnames(data_to_append), ]
            if(class(col_data_true) %in% "character") col_data_true <- data.frame(id=col_data_true)
          	col_data_true <- col_data_true[order(col_data_true$id), ]
            if(class(col_data_true) %in% "character") col_data_true <- data.frame(id=col_data_true)
          	cluster_data_true <- colnames(data_to_append)[grep("eigen*", colnames(data_to_append))]
          	cluster_data_true <- data.frame(id=cluster_data_true)
          	rownames(cluster_data_true) <- cluster_data_true$id
          
          	col_data_true <- plyr::rbind.fill(col_data_true, cluster_data_true)
          	col_data_true <- col_data_true[order(col_data_true$id), ]
          	col_info <- col_info[order(col_info$id), ]
          	col_data_true <- cbind.data.frame(col_data_true, col_info$modules)
            colnames(col_data_true)[length(colnames(col_data_true))] <- "modules"
          

            name <- paste(which_data, "with_eigenmetabs_1", sep="_")
            while(name %in% names(object@list_of_data)) {
                index <- strsplit(name, split="_") %>% lapply(function(x) return(x[length(x)])) %>% unlist() %>% as.numeric()
                index <- index + 1
                name <- paste(which_data, "with_eigenmetabs", index, sep="_")
            }

          	if(append) {
          		out <- add_dataset(object, data=data_to_append, col_data=col_data_true, 
                    row_data=object@list_of_row_data[[which_data]], name=name)
          	} else {
          		out <- get_make_analyser_object(data=data_to_append, col_data=col_data_true, 
          				row_data=object@list_of_row_data[[whichd_ata]], 
          				annotations_index=list(), name=paste(which_data, "with_eigenmetabs", sep="_"), results=list())
          	}
            out <- add_function_info(object=out, function_name="mod_trans_eigendata", 
                params=list(which_data=which_data, append=append, results_index=results_index, cols_for_meta=cols_for_meta, name=name, ...))
          	return(out)
	})