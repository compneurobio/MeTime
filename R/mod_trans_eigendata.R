#' Function to add the clusters obtained from wgcna
#' @description A function to add cluster assiginment to the col_data from WGCNA
#' @param object An S4 object of class metime_analyser
#' @param which_data character to define which dataset is to be used
#' @param append logical if set to true adds the new data to the object used else creates new object
#' @param clusters logical if set to true will add already existing cluster info otherwise creates new
#' @return metime_analyser object with new dataset with eigendata of the metabolites
#' @export
setGeneric("mod_trans_eigendata", function(object, which_data, append, clusters, ...) standardGeneric("mod_trans_eigendata"))
setMethod("mod_trans_eigendata", "metime_analyser", function(object, which_data, append, clusters, ...) {
			stopifnot(which_data %in% names(object@list_of_data))
      if(!clusters) {
          object <- add_clusters_wgcna(object=object, which_data=which_data, ...)
      }
			cluster_info <- data.frame(id=object@list_of_col_data[[which_data]]$id,
										cluster=object@list_of_col_data[[which_data]][ ,grep("wgcna_clusters", colnames(object@list_of_col_data[[which_data]]))])
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

            col_data <- get_metadata_for_columns(object=object, which_data=which_data, columns=list(c("id", "sub_pathway")), 
                 names=c("id", "pathway"), index_of_names=rep("id", each=length(which_data)))
          	col_data_true <- col_data[col_data$id %in% colnames(data_to_append), ]
          	col_data_true <- col_data_true[order(col_data_true$id), ]
          	col_data_true <- col_data_true[,1:2]
          	colnames(col_data_true)[2] <- "pathway"
          	cluster_data_true <- colnames(data_to_append)[grep("eigen*", colnames(data_to_append))]
          	cluster_data_true <- data.frame(id=cluster_data_true, pathway=cluster_data_true)
          	rownames(cluster_data_true) <- cluster_data_true$id
          
          	col_data_true <- rbind.data.frame(col_data_true, cluster_data_true)
          	col_data_true <- col_data_true[order(col_data_true$id), ]
          	col_info <- col_info[order(col_info$id), ]
          	col_data_true <- cbind.data.frame(col_data_true, col_info$modules)
          	colnames(col_data_true)[3] <- c("modules")
          
          	if(append) {
          		out <- get_append_analyser_object(object, data=final_data_true, col_data=col_data_true, 
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