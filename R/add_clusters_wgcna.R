#' Function to add the clusters obtained from wgcna
#' @description A function to add cluster assiginment to the col_data from WGCNA
#' @param object An S4 object of class metime_analyser
#' @param which_data character to define which dataset is to be used
#' @param baseline character to define the timepoint to be used as baseline to predict clusters
#' @param ... other parameters for cutreeDynamic function such minClusterSize, pamDendroRespect etc
#' @return metime_analyser object with updated column info about the clustersize
#' @export
setGeneric("add_clusters_wgcna", function(object, which_data, baseline, ...) standardGeneric("add_clusters_wgcna"))
setMethod("add_clusters_wgcna", "metime_analyser", function(object, which_data, baseline, ...) {
		data <- object@list_of_data[[which_data]]
        col_data <- object@list_of_col_data[[which_data]]
        data <- data[grep(baseline, rownames(data)), ]
        lnames <- names(data)
        # Choose a set of soft-thresholding powers
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        # Call the network topology analysis function
        sft = WGCNA::pickSoftThreshold(data, powerVector = powers, verbose = 5)
        res <- sft$fitIndices[,1][which(-sign(sft$fitIndices[,3])*sft$fitIndices[,2] > 0.88)] 
        softPower <- res[1]
        adjacency <- WGCNA::adjacency(data, power=softPower)
        # Turn adjacency into topological overlap
        TOM = WGCNA::TOMsimilarity(adjacency)
        dissTOM = 1-TOM
        # Call the hierarchical clustering function
        geneTree = hclust(as.dist(dissTOM), method = "average")  
        dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM, ...)
        modules <- cbind.data.frame(lnames, dynamicMods)
        modules <- modules[order(modules$lnames), ]
        object@list_of_col_data[[which_data]]$wgcna_clusters <- modules[ ,2]
        colnames(object@list_of_col_data[[which_data]])[colnames(object@list_of_col_data[[which_data]]) %in% "wgcna_clusters"] <- paste("wgcna_clusters", ..., sep="_")
        object <- add_function_info(object=object, 
                        function_name="add_clusters_wgcna", 
                        params=list(which_data=which_data, baseline=baseline, ...))
        return(object)
	})