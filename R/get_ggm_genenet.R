
#' Function to calculate a dynamic GeneNet GGM from a longitudnal data matrix
#' @description calculates GGM on longitudnal data matrix and returns a dataframe with edges, 
#'   partial correlation and associated p-values
#' @param data data matrix in a longitudnal format
#' @param threshold type of multiple hypothesis correction. Available are Bonferoni("bonferroni"), 
#'   Benjamini-Hochberg("FDR") and independent tests method("li", also see Li et al ....)
#' @param all Logical to get all edges without any cutoff.
#' @param ... additional arguments for ggm.estimate.pcor()
#' @return a dataframe with edges, partial correlation and associated p-values 
#' @export
get_ggm_genenet <- function(data, threshold=c("bonferroni", "FDR", "li"), all, ...) {
  # check if longitudinal
  if(!longitudinal::is.longitudinal(data)) stop("data is not a longitudinal object") 
  
  met.ggm <- GeneNet::ggm.estimate.pcor(data, method="dynamic", ...) # retrieve GGM
  met.ggm.edges <- GeneNet::network.test.edges(met.ggm, plot=F) # calculate edge statistics
  
  #define thresholds
  p.thresh <- 0.05/((ncol(met.ggm))*(ncol(met.ggm))/2) 
  fdr.thresh <- 0.05
  #Check all or threshold
  if(all) {
  	met.ggm.edges.filtered <- met.ggm.edges
  } else {
  	# cut at threshold
  	if(threshold=="FDR") {
      met.ggm.edges.filtered <- met.ggm.edges[which(met.ggm.edges$qval < 0.05),]
  	} else if(threshold=="li"){
      data <- data %>% as.matrix() %>% .[,] %>% as.data.frame()  
      cordat <- cor(data)
      eigenvals <- eigen(cordat)$values
      li.thresh <- sum( as.numeric(eigenvals >= 1) + (eigenvals - floor(eigenvals)) )
      met.ggm.edges.filtered <- met.ggm.edges[which(met.ggm.edges$pval < 0.05/li.thresh),]
 	 	} else if(threshold=="bonferroni"){
      met.ggm.edges.filtered <- met.ggm.edges[which(met.ggm.edges$pval < p.thresh),]
  	}
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
  for(i in 1:nrow(met.ggm.edges.filtered)) {
    cor.nodes <- cor.test(data[,met.ggm.edges.filtered$node1[i]],data[,met.ggm.edges.filtered$node2[i]])
    # Print and store those that do not make it
    if(cor.nodes$p.value > p.thresh){
      cat(met.ggm.edges.filtered$node1[i]," : ", met.ggm.edges.filtered$node2[i], " -> pcor=", met.ggm.edges.filtered$pcor[i],"(P=",met.ggm.edges.filtered$pval[i],"), cor=", cor.nodes$estimate, "(P=", cor.nodes$p.value,")\n")
      edge2rem <- c(edge2rem, i)
    }
  }
  
  # Remove edges without significant pairwise correlations
  out <- met.ggm.edges.filtered[-edge2rem,]
  return(out)
}

