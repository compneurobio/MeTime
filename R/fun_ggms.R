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
adni_ggm_calc_ggm_dynamic <- function(data, threshhold=c("bonferoni", "FDR", "li")){
  # check if longitudinal
  if(!longitudinal::is.longitudinal(data)) stop("data is not a longitudinal object")
  
  met.ggm <- GeneNet::ggm.estimate.pcor(data, method="dynamic")     # retrieve GGM
  met.ggm.edges <- GeneNet::network.test.edges(met.ggm, plot=F)           # calculate edge statistics
  
  #define thresholds
  p.thresh <- 0.05/((ncol(met.ggm))*(ncol(met.ggm))/2) 
  fdr.thresh <- 0.05
  li.thresh <- data %>% as.matrix() %>% 
    .[,] %>% as.data.frame() %>% 
    adni_independent_tests()
    
  # cut at threshold
  if(threshhold=="FDR"){
    met.ggm.edges.filtered <- met.ggm.edges[which(met.ggm.edges$qval < 0.05),]
  }
  else if(threshhold=="li"){
    met.ggm.edges.filtered <- met.ggm.edges[which(met.ggm.edges$pval < 0.05/li.thresh),]
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

  
# data<- database$lipid$data %>% 
#   adni_add_index() %>%  # add index
#   adni_filter_full_tp(tp=c(0,12,24)) %>%  # filter all subjects measured at all three time points
#   adni_ggm_convert_longitudinal() %>% # convert to longitudinal object
#   adni_ggm_calc_ggm_dynamic(threshhold="li")

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

}

#' Function to plot data from network and object after calculating a certain ggm
#' @description A function to plot ggms of different kinds namely visNetwork, cytoscape and networkD3
#' @param type Type of visualization options = c("visNetwork", "Cytoscape", "networkD3")
#' @param network dataframe with information of the network
#' @param metadata dataframe with two columns 1) metabs and 2) Group they belong to
#' @return network plot based on the type chosen by the user
ggm_visualizer <- function(network, type, metadata) {
      #sanitychecks
      stopifnot(colnames(metadata) %in% c("name","group","class"))
      stopifnot(type %in% c("visNetwork", "Cytoscape", "networkD3"))
      if(type %in% "networkD3") {

      } else if(type %in% "visNetwork") {

      } else {

      }
}

#' An automated function to caluclate temporal network with lagged model
#' @description calculates temporal networks for each dataset with a lagged model as used in graphical VAR
#' @param object object obtained after applying prep_data_for_ggms() on S4 object of class metab_analyser
#' @param lag which lagged model to use. 1 means one-lagged model, similary 2,3,..etc
#' @param timepoints timepoints of interest that are to be used to build networks
#' @return temporal network data with edgelist and regression values
#' @export
automated_temporal_network <- function(object, lag, timepoints) {
    
}

#USAGE
#prep data like always
#object <- prep_data_for_ggms(object)
#network <- automated_ggm_genenet(object, which_data, threshold, timepoints)
#network <- automated_ggm_mlp(object, which_data, rho, timepoints)
#network_plot <- ggm_visualizer(object, metadata, type)