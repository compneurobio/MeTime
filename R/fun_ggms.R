# script to calculate ggms




#' Function to convert data into a longitudnal format for GeneNet ggms
#' @description converts a dataset with full data into a longitudnal version
#' @param data data matrix with metabolite concentrations
#' @return data matrix converted into a longitudnal format
#' @export
adni_ggm_convert_longitudinal <- function(data){
  # get all variable names
  vars <- data %>% 
    adni_rm_index() %>% 
    names()
  # get full data
  data <- data %>% 
    adni_add_index() %>% 
    dplyr::arrange(timepoint, rid) 
  n_subject = unique(data$rid) %>% length() %>% as.numeric()
  name_tp = unique(data$timepoint) %>% as.numeric()
  
  out <- longitudinal::as.longitudinal(
    x = as.matrix(data[,vars]),
    repeats = n_subject,
    time = name_tp
  )
  return(out)
}

#check this now:
#object$list_of_data <- lapply()
#object$list_of_data <- lapply(object$list_of_data, adni_ggm_convert_longitudnal) 


#' Function to calculate a dynamic GeneNet GGM
#' @description calculates GGM on longitudnal data matrix and returns a dataframe with edges, partial correlation and associated p-values
#' @param data data matrix in a longitudnal format(see adni_convert_ggm_longitudnal)
#' @param threshold type of multiple hypothesis correction. Available are Bonferoni("bonferoni"), Benjamini-Hochberg("FDR") and independent tests method("li", also see Li et al ....)
#' @return a dataframe with edges, partial correlation and associated p-values 
#' @export
adni_ggm_calc_ggm_dynamic <- function(data, threshhold=c("bonferoni", "FDR", "li")){
  # check if longitudinal
  if(!longitudinal::is.longitudinal(data)) stop("data is not a longitudinal object")
  
  met.ggm <- GeneNet::ggm.estimate.pcor(data, method="dynamic")     # retrieve GGM
  met.ggm.edges <- GeneNet::network.test.edges(met.ggm, plot=F)           # calculate edge statistics
  
  #define thresholds
  p.thresh <- 0.05/((ncol(met.ggm))*(ncol(met.ggm))/2) 
  fdr.thresh <- 0.05
    
  # cut at threshold
  if(threshhold=="FDR"){
    met.ggm.edges.filtered <- met.ggm.edges[which(met.ggm.edges$qval < 0.05),]
  }
  else if(threshhold=="li"){
    data <- data %>% as.matrix() %>% .[,] %>% as.data.frame()  
    cordat <- cor(data)
    eigenvals <- eigen(cordat)$values
    li.thresh <- sum( as.numeric(eigenvals >= 1) + (eigenvals - floor(eigenvals)) )
    met.ggm.edges.filtered <- met.ggm.edges[which(met.ggm.edges$pval < 0.05/li.thresh),]
}

  }
  else if(threshhold=="bonferoni"){
    met.ggm.edges.filtered <- met.ggm.edges[which(met.ggm.edges$pval < p.thresh),]
  }
  
  ## Reinsert node (metabolite) names
  node1list <- NULL
  node2list <- NULL
  for(i in 1:nrow(met.ggm.edges.filtered)){
    node1list <- c(node1list, colnames(met.ggm)[met.ggm.edges.filtered$node1[i]])
    node2list <- c(node2list, colnames(met.ggm)[met.ggm.edges.filtered$node2[i]])
  }
  met.ggm.edges.filtered$node1 <- node1list
  met.ggm.edges.filtered$node2 <- node2list
  ## Filter edges for significant partial correlations that are also significant pairwise correlations
  edge2rem <- NULL
  for(i in 1:nrow(met.ggm.edges.filtered)){
    cor.nodes <- cor.test(data[,met.ggm.edges.filtered$node1[i]],data[,met.ggm.edges.filtered$node2[i]])
    # Print and store those that do not make it
    if(cor.nodes$p.value > p.thresh){
      cat(met.ggm.edges.filtered$node1[i]," : ", met.ggm.edges.filtered$node2[i], " -> pcor=", met.ggm.edges.filtered$pcor[i],"(P=",met.ggm.edges.filtered$pval[i],"), cor=", cor.nodes$estimate, "(P=", cor.nodes$p.value,")\n")
      edge2rem <- c(edge2rem, i)
    }
  }
  
  # Remove edges without significant pairwise correlations
  met.ggm.edges.filtered <- met.ggm.edges.filtered[-edge2rem,]
  return(met.ggm.edges.filtered)
}


#' An automated fucntion to calculate GGM from genenet
#' @description automated funtion that can be applied on s3 object obtained after prep_data_for_ggms() to obtain geneNet network along with threshold used
#' @param object object obtained after applying prep_data_for_ggms() on S4 object of class metab_analyse with type="single"
#' @param which_data a character or a character vector naming the datasets of interest
#' @param threshold type of threshold to be used for extracting significant edges
#' @param timepoints timepoints of interest that are to be used to build networks
#' @return Network data with edgelist, partial correlation values and associated p-values and corrected p-values 
automated_ggm_genenet <- function(object, which_data, threshhold, timepoints) {
    #sanity checks
    stopifnot(timepoints %in% c("t0","t12","t24"))
    stopifnot(threshhold %in% c("li","bonferoni","FDR"))
    #Extracting data that is needed
    col_data_names <- gsub("_data", "_col_data", which_data)
    row_data_names <- gsub("_data", "_row_data", which_data)
    object$list_of_data <- object$list_of_data[names(object$list_of_data) %in% which_data]
    object$list_of_col_data <- object$list_of_col_data[names(object$list_of_col_data) %in% col_data_names]
    object$list_of_row_data <- object$list_of_row_data[names(object$list_of_row_data) %in% row_data_names]
    timepoints <- as.character(unlist(lapply(strsplit(timepoints), split="t", fixed=TRUE), function(x) return(x[2])))
    #converting prepped data to full ggm network
    if(length(which_data > 1)) {
          data <- do.call(cbind, object$list_of_data)
          colnames(data) <- unlist(lapply(strsplit(colnames(data), split="_data."), function(x) return(x[2])))
          network <- data %>%
                    adni_add_index() %>%
                    adni_filter_full_tp(tp = timepoints) %>%
                    adni_ggm_calc_ggm_dynamic(threshhold = threshhold) 
      } else {
          network <- object$list_of_data[[1]]
          network <- network %>%
                    adni_add_index() %>%
                    adni_filter_full_tp(tp = timepoints) %>%
                    adni_ggm_calc_ggm_dynamic(threshhold = threshhold)
         
      }
    return(network)
    #Add visualization function here
}

#' An automated fucntion to calculate GGM from genenet
#' @description automated funtion that can be applied on s3 object obtained after prep_data_for_ggms() to obtain geneNet network along with threshold used. 
#' This function is not applicable for singular datasets
#' @param object object obtained after applying prep_data_for_ggms() on S4 object of class metab_analyse with type="multi" not applicable to apply on singular sets 
#' @param which_data a character or a character vector naming the datasets of interest
#' @param rho tuning parameter for regression
#' @param nfolds check
#' @param timepoints timepoints of interest that are to be used to build networks
#' @return Network data with edgelist, partial correlation values and associated p-values and corrected p-values 
automated_ggm_mlp <- function(object, which_data, rho, nfolds, timepoints) {
    #Prepping data for performing regression
    object$list_of_data <- object$list_of_data[names(object$list_of_data) %in% which_data]
    object$list_of_data <- lapply(object$list_of_data, function(x) {
          x <- x[names(x) %in% timepoints]
          return(x) 
    })
    if(length(which_data) > 1) {
        data_list_time <- lapply(object$list_of_data, function(x) {
              x <- lapply(x, function(r) return(r[sort(rownames(r)),]))
          })
        final_data_list <- list()
        count <- 1
        for(i in 1:length(data_list_time[[1]])) {
          foo <- lapply(data_list_time, function(x) return(x[i]))
          final_data_list[[count]] <- do.call(cbind, foo)
        } 
        names(final_data_list) <- names(data_list_time[[1]])  
    } else {
        final_data_list <- object$list_of_data
    }
}

#' Function to plot data from network and object after calculating a certain ggm
#' @description A function to plot ggms of different kinds namely visNetwork and cytoscape. Make sure that cytoscape is running before running this function
#' @param type Type of visualization options = c("visNetwork", "Cytoscape")
#' @param network dataframe with information of the network colnames should be in this format node1, node2, corr_val and other columns 
#' @param metadata dataframe with three columns and another optional column 1) metabs("name") and 2) Group("group") they belong to 3) and the class("class") they belong to
#' see get_metadata_for_plotting() for more information. The optional column is for colors
#' @param type_of_data character defining the type of network applied to a single network or a multi dataset network. Available options = c("single","multi")
#' @return network plot based on the type chosen by the user
#Add vector option for colors and timepoint fold change network
#check how many metabs are matching back to metadata 
#add layout parameters
ggm_visualizer <- function(network, type_of_plot, metadata, main, type_of_data, timepoints_fold, layout_by) {
      
      #sanitychecks
      stopifnot(colnames(metadata) %in% c("name","group","class"))
      stopifnot(type %in% c("visNetwork", "Cytoscape"))
      
      #Graph building based on the type of graph requested 
      #Graph for visNetwork based network
      if(type_of_plot %in% "visNetwork") {
          #Preparing data
          #Getting Node list
          nodes <- unique(c(network$node1, network$node2))
          node_list <- data.frame(id=1:length(nodes), label=nodes, group=as.character(1:length(nodes)))
          for(i in 1:length(node_list$label)) {
            g <- metadata[as.character(metadata$name) %in% as.character(node_list$label[i]), "group"]
            node_list$group[i] <- g
          }
          
          #Getting edge list
          edge_list <- data.frame(from=1:length(network$node1), to=1:length(network$node2))
          for(i in 1:length(network$node1)) {
            edge_list$from[i] <- node_list[as.character(node_list$label) %in% as.character(data$node1[i]), "id"]
            edge_list$to[i] <- node_list[as.character(node_list$label) %in% as.character(data$node2[i]), "id"]
          }
          dashes <- ifelse(data$pcor_val > 0, FALSE, TRUE)
          edge_list$dashes <- dashes
          edge_list$values <- data$pcor_val
          
          #Choosing colors and shapes for visualization
          if(is.null(metadata$colors)) {
              shapes <- c("square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond")
              colors_for_nodes <- metadata$group
              colors_code <- cbind(unique(node_list$group), get_palette(length(unique(node_list$group))))
              colnames(colors_code) <- c("group", "color")
              for(i in 1:length(colors_for_nodes)) {
                  colors_for_nodes[i] <- colors_code[colors_for_nodes %in% colors_code$group, "color"]
              }
          } else {
              shapes <- c("square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond")
              color_for_nodes <- metadata$colors
          }
          
          ledges <- data.frame(color = c("#920000","#0072b2"), label = c("negative", "positive"), dashes =c(TRUE, FALSE))
          
          #initiating the graph based on type of data - multi
          if(type_of_data %in% "multi") { 
              classes <- unique(metadata$class)
              groups <-  unique(metadata$group)
              graph <- visNetwork(nodes=node_list, edges=edge_list, main=main) %>%
                    visIgraphLayout(layout=layout_by, physics = F, smooth = F) %>%
                    visPhysics(stabilization = FALSE)
              for(i in 1:length(classes)) {
                 graph <- visGroups(graph=graph, groupname = classes[i], shape = shapes[i])
              }
              graph <- visEdges(graph=graph, addEdges = ledges, useGroups = T) %>% 
                        visNodes(borderWidth = 3, color=list(background=colors_for_nodes))
                        visOptions(highlightNearest = list(enabled=T, hover=T), nodesIdSelection = T, selectedBy = "group")%>%
                        visInteraction(navigationButtons = T) %>%
                        visConfigure(enabled=T)

          #initiating the graph based on type of data - multi
          } else if(type_of_data %in% "single") {
              graph <- visNetwork(nodes=node_list, edges=edge_list, main=main) %>%
                    visIgraphLayout(layout=layout_by, physics = F, smooth = F) %>%
                    visPhysics(stabilization = FALSE) %>%
                    visEdges(addEdges = ledges, useGroups = T) %>% 
                    visNodes(borderWidth = 3, color=list(background=colors_for_nodes))
                    visOptions(highlightNearest = list(enabled=T, hover=T), nodesIdSelection = T, selectedBy = "group")%>%
                    visInteraction(navigationButtons = T) %>%
                    visConfigure(enabled=T)

          } else {
               stop("Check the input for type_of_data. The type chosen is not available")
          }
          
      #Graph for cytoscape based network
      } else if(type_of_plot %in% "cytoscape") {
          #prepping data i.e. getting edge and node tables
          edge <- network[,1:3]
          colnames(edge) <- c("source", "target", "weight")
          edge$interaction <- ifelse(edge$weight > 0, "positive", "negative")
          edge$weight <- abs(edge$weight) 
          graph <- igraph::simplify(igraph::graph.data.frame(edge, directed=FALSE))
          node <- data.frame(id=igraph::V(graph)$name)
          node$group <- rep("g", each=length(node$id))
          for(i in 1:length(node$id)) {
            g <- metadata[as.character(metadata$name) %in% as.character(node$id[i]), "group"]
            node$group[i] <- g
          }
          
          #Quantifying betweenness for all nodes and adding those to node dataframe
          bet_all <- igraph::betweenness(graph, v = igraph::V(graph), directed = FALSE) / (((igraph::vcount(graph) - 1) * (igraph::vcount(graph)-2)) / 2)
          bet_all_norm <- (bet_all - min(bet_all))/(max(bet_all) - min(bet_all))
          node <- as.data.frame(cbind(node, score=100*betAll.norm))

          #Changing factors in edge table to characters.  
          cols_with_factors <- sapply(edge, is.factor)
          edge[cols_with_factors]<- lapply(edge[cols_with_factors], as.character)

          #initiating the graph based on type of data - single
          if(type_of_data %in% "single") {
              graph <- createNetworkFromDataFrames(node, edge, title=main, collection="DataFrame_type") %>%
                        setNodeColorMapping("group", unique(node$group), get_palette(length(unique(node$group)), mapping.type="d") %>%
                        setEdgeLineWidthMapping('weight', c(max(edge$weight), min(edge$weight)), c(10*max(edge$weight), 10*min(edge$weight))) %>%
                        setEdgeLineStyleMapping("interaction", c("positive", "negative"), c("SOLID","LONG_DASH")) %>%
                        setNodeCustomLinearGradient(c("#DDDDDD", "#888888"),anchors = c(10*min(node$score), 10*max(node$score)))

          #initiating the graph based on type of data - multi
          } else if(type_of_data %in% "multi") {
              #Adding shapes to showcase different kinds of data
              shapes <- getNodeShapes()
              node$class <- rep("g", each=length(node$id))
              for(i in 1:length(node$id)) {
                  g <- metadata[as.character(metadata$name) %in% as.character(node$id[i]), "class"]
                  node$class[i] <- g
              }
              #Graph
              graph <- createNetworkFromDataFrames(node,edge, title=main, collection="DataFrame_type") %>%
                        setNodeColorMapping("group", unique(node$group), get_palette(length(unique(node$group)), mapping.type="d"), mapping.type="d") %>%
                        setNodeShapeMapping("class", unique(node$class), shapes[length(unique(node$class))]) %>%
                        setEdgeLineWidthMapping('weight', c(max(edge$weight), min(edge$weight)), c(10*max(edge$weight), 10*min(edge$weight))) %>%
                        setEdgeLineStyleMapping("interaction", c("positive", "negative"), c("SOLID","LONG_DASH")) %>%
                        setNodeCustomLinearGradient(c("#DDDDDD", "#888888"),anchors = c(10*min(node$score), 10*max(node$score))) %>%
          } else {
               stop("Check the input for type_of_data. The type chosen is not available")
          }

      } else {
         stop("Check the input for type_of_plot. The type chosen is not available")
      }
      return(graph)
}

#' An automated function to caluclate temporal network with lagged model
#' @description calculates temporal networks for each dataset with a lagged model as used in graphical VAR
#' @param object object obtained after applying prep_data_for_ggms() on S4 object of class metab_analyser
#' @param lag which lagged model to use. 1 means one-lagged model, similary 2,3,..etc
#' @param which_data dataset or datasets to be used
#' @param timepoints timepoints of interest that are to be used to build networks(in the order of measurement)
#' @param rho parameter for regression coefficient
#' @return temporal network data with edgelist and regression values
#' @export
automated_temporal_network <- function(object, lag, rho, timepoints, which_data) {
    object$list_of_data <- object$list_of_data[names(object$list_of_data) %in% which_data]
    object$list_of_data <- lapply(object$list_of_data, function(x) {
          x <- x[names(x) %in% timepoints]
          return(x) 
    })
    number_of_models <- length(timepoints) - lag 
    count <- 1
    model_seqs <- list()
    while(lag + count <= length(timepoints)) {
      model_seqs[[count]] <- timepoints[seq(from=count, to=lag+count, by=1)]
      count <- count + 1
    }
    if(length(which_data) > 1) {
        data_list_time <- lapply(object$list_of_data, function(x) {
              x <- lapply(x, function(r) return(r[sort(rownames(r)),]))
          })
        final_data_list <- list()
        count <- 1
        for(i in 1:length(data_list_time[[1]])) {
          foo <- lapply(data_list_time, function(x) return(x[i]))
          final_data_list[[count]] <- do.call(cbind, foo)
        } 
        names(final_data_list) <- names(data_list_time[[1]])  
    } else {
        final_data_list <- object$list_of_data
    }
    
    models <- list()
    
}

#' Function to perform regression on list of matrics either divided based on time(temporal net) or data type(Multibipartite Lasso)
#' @param list_of_mats list of matrices that are divided based on platform or timepoints
#' @param alpha parameter for glmnet alpha=1 represents ridge regression and alpha=0 represents lasso regression and anything in between results in a mixed penalty regression
#' @param nfolds nfolds parameter for glmnet
#' @return a list with different combinations used and each combination is a nested list with regression result data
#' @export
get_betas <- function(list_of_mats, 
           alpha,
           nfolds 
           ) {
  #creating a list to store the data from glmnet
  #code exactly similar to the usual MLP 
  list_with_combos <- list() # list to store the regression information for each metabolte
  count <- 1
  for(i in 1:length(list_of_mats)) {
    for(j in 1:length(list_of_mats)) {
      if(i != j) {
        x_mat <- as.matrix(list_of_mats[[i]])
        y_mat <- as.matrix(list_of_mats[[j]])
        fit_list <- list()
        for(k in 1:ncol(y_mat)) {
          fit_list[[k]] <- cv.glmnet(x=x_mat, y=y_mat[,k], alpha=alpha, nfolds=nfolds)
          names(fit_list)[k] <- colnames(y_mat)[k]
        }
        list_with_combos[[count]] <- fit_list
        names(list_with_combos)[count] <- paste(names(list_of_mats)[j], names(list_of_mats)[i], sep="-")
        count <- count+1
      } else {
        fit_list <- list()
        mat <- as.matrix(list_of_mats[[i]])
        for(k in 1:ncol(mat)) {
          y <- as.matrix(mat[,k])
          x <- as.matrix(mat[,-k])
          fit_list[[k]] <- cv.glmnet(x=x, y=y, alpha=alpha, nfolds=nfolds)
          names(fit_list)[k] <- colnames(mat)[k]
        } 
        list_with_combos[[count]] <- fit_list
        names(list_with_combos)[count] <- names(list_of_mats)[i]
        count <- count +1
      }

    }
  }
  return(list_with_combos)
}



#USAGE
#prep data like always
#object <- prep_data_for_ggms(object, which_type, mlp_or_temp)
#network <- automated_ggm_genenet(object, which_data, threshold, timepoints)
#network <- automated_ggm_mlp(object, which_data, rho, timepoints)
#network_plot <- ggm_visualizer(network, metadata, type)



######ROUFNFFNFNkdsgknognpsfvlma##sdngksnflm#vS

#m <- lapply(k, function(a) {
#    a <- lapply(a, function(r) {
#        return(r[sort(rownames(r)),])
#      })
#  })