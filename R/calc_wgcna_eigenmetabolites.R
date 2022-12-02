
#' Function to extract eigenmetabolites and modules from WGCNA
#' @description Function to calculate modules and eigenmetabolites
#' @param object An s4 object of class metime_analyser
#' @param which_data Dataset to be used for the analysis
#' @param baseline baseline timepoint
#' @param minClusterSize minimum number of metabolites in a particular cluster
#' @return plotter object with dendogram results
#' @export
setGeneric("calc_wgcna_eigenmetabolites", function(object, which_data, baseline, minClusterSize, ...) standardGeneric("calc_wgcna_eigenmetabolites"))
setMethod("calc_wgcna_eigenmetabolites", "metime_analyser", function(object, which_data, baseline, ...) {
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
        if(is.null(minClusterSize)) {
          dynamicMods_w_respect = dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM,
                    deepSplit = 3, pamRespectsDendro = TRUE)
        } else {
          dynamicMods_w_respect = dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM,
                    deepSplit = 3, pamRespectsDendro = TRUE, minClusterSize=minClusterSize)
        }
    
        #module data
        modules_w_respect <- cbind.data.frame(lnames, dynamicMods_w_respect)
        modules_w_respect <- modules_w_respect[order(modules_w_respect$lnames), ]
        object@list_of_col_data[[which_data]]$cluster_wgcna_true <- modules_w_respect[ ,2]

        #splitting data into modules
        data_list_with_respect <- split(modules_w_respect, f = modules_w_respect[ ,2])

        #making new data with eigenmetabolite calculation for each cluster
        newdata <- object@list_of_data[[which_data]]
        final_data_true <- lapply(names(data_list_with_respect), function(x) {
                  modules <- data_list_with_respect[[x]]
                  if(x %in% "0") {
                      dummy <- newdata[ ,modules[,1]]
                  } else {
                      dummy <- newdata[ ,modules[ ,1]]
                      pca <- prcomp(dummy, scale.=T, center=T)
                      dummy <- as.data.frame(pca$x[,1])
                      colnames(dummy) <- c(paste("eigenmetabolite_cluster_", x, sep=""))
                  }
                  return(dummy)
          }) %>% do.call(what=cbind.data.frame)

          col_data <- get_metadata_for_columns(object=object, which_data=which_data, columns=list(c("id", "sub_pathway")), 
                 names=c("id", "pathway"), index_of_names=rep("id", each=length(which_data)))
          col_data_true <- col_data[col_data$id %in% colnames(final_data_true), ]
          col_data_true <- col_data_true[order(col_data_true$id), ]
          col_data_true <- col_data_true[,1:2]
          colnames(col_data_true)[2] <- "pathway"
          truedata <- data.frame(table(col_data_true[ ,2]))
          truedata <- cbind.data.frame(truedata, rep("True", each=length(truedata[,1])))
          colnames(truedata) <- c("group", "count", "type")
          
          barplot <- ggplot(truedata, aes(x=group, y=count, fill=type)) + geom_bar(stat="identity", position=position_dodge()) + coord_flip()

          cluster_data_true <- colnames(final_data_true)[grep("eigen*", colnames(final_data_true))]
          cluster_data_true <- data.frame(id=cluster_data_true, pathway=cluster_data_true)
          rownames(cluster_data_true) <- cluster_data_true$id
          
          col_data_true <- rbind.data.frame(col_data_true, cluster_data_true)
          
          object <- get_append_analyser_object(object, data=final_data_true, col_data=col_data_true, 
                    row_data=object@list_of_row_data[[which_data]], name=paste(which_data, "with_eigenmetabs", sep="_"))
  
        #table(dynamicMods)
        #dynamicColors <- labels2colors(dynamicMods)
        #plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
        #        dendroLabels = FALSE, hang = 0.03,
        #        addGuide = TRUE, guideHang = 0.05,
        #        main = "Metabolites dendrogram and module colors")  

        return(list(object=object, plot=barplot))
  })


