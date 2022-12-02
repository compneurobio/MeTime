
#' An automated fucntion to calculate GGM from genenet crosssectional version
#' @description automated funtion that can be applied on metime_analyser object to obtain geneNet network along with threshold used
#' @param object S4 object of cĺass metime_analyser
#' @param which_data a character or a character vector naming the datasets of interest
#' @param threshold type of threshold to be used for extracting significant edges. 
#'      allowed inputs are "li", "FDR", "bonferroni"
#' @param timepoints timepoints of interest that are to be used to build networks(as per timepoints in rows)
#' @param all Logical to extract all edges without any pval correction
#' @param cols_for_meta list of character vector for extracting metadata of metabolites for plotting
#' @param covariates covariates to be used for this analysis
#' @param ... additional arguments for GeneNet
#' @return Network data as a plotter object
#' @export
setGeneric("calc_ggm_genenet_crosssectional", function(object, which_data, threshold, timepoint, all, cols_for_meta, covariates, ...) standardGeneric("calc_ggm_genenet_crosssectional"))
setMethod("calc_ggm_genenet_crosssectional", "metime_analyser", function(object, which_data, threshold, timepoint, all, cols_for_meta, covariates, ...) {
        if(length(which_data) > 1) object <- mod_extract_common_samples(object)
        stopifnot(threshold %in% c("li","bonferroni","FDR", NULL))
        object <- mod_split_acc_to_time(object)
        list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
        list_of_data <- lapply(list_of_data, function(x) return(x[names(x) %in% timepoint]))
        list_of_data <- lapply(list_of_data, function(x) return(x[[1]][order(rownames(x[[1]])),]))
        list_of_data <- unname(list_of_data)
        data <- as.data.frame(do.call(cbind, list_of_data))
        rm_col = intersect(colnames(data), c("adni_id","RID","rid","timepoint","tp","subject", "id"))
        data <- data %>% select(-c(all_of(rm_col)))
        data <- na.omit(data)
        this_mat <- as.matrix(apply(data, 2, as.numeric))
        pcor_mat <- GeneNet::ggm.estimate.pcor(as.matrix(this_mat), method = "dynamic", verbose = F)
        # compute p-values of edges
        pval_mat <- GeneNet::network.test.edges(pcor_mat, plot = F, verbose = F)
        # p-value correction of edges
        pval_mat$p.adj.bh <- p.adjust(pval_mat$pval, method="BH")
        pval_mat$p.adj.bon <- p.adjust(pval_mat$pval, method="bonferroni")    
        ggm_thresh <- 0.05
        # extract edge list
        tmp <- pcor_mat %>% igraph::graph_from_adjacency_matrix(mode='undirected', weighted = T) %>% igraph::simplify()
        ggm_edges <- cbind.data.frame(get.edgelist(tmp), edge_attr(tmp)$weight)
        names(ggm_edges) <- c("node1", "node2", "pcor")
        if(all) {
           ggm_data <- ggm_edges
        } else {
          if(threshold %in% "bonferroni") {
            #filter edges based on p-values - bonferroni
            ggm_data <-  ggm_edges %>% filter(abs(pcor)>=min(abs(pval_mat$pcor[pval_mat$p.adj.bon<=ggm_thresh]))) 
          } else if(threshold %in% "FDR"){
           # filter edges based on p-values - BH 
            ggm_data <-  ggm_edges %>% filter(abs(pcor)>=min(abs(pval_mat$pcor[pval_mat$p.adj.bh<=ggm_thresh])))
          } else if(threshold %in% "li") {
            data <- this_mat %>% as.matrix() %>% .[,] %>% as.data.frame()  
            cordat <- cor(data)
            eigenvals <- eigen(cordat)$values
            li.thresh <- sum( as.numeric(eigenvals >= 1) + (eigenvals - floor(eigenvals)) )
            ggm_data <- ggm_edges %>% filter(abs(pcor)>=min(abs(pval_mat$pcor[pval_mat$pval<=0.05/li.thresh])))
          }  
        }
    network <- ggm_data
    network <- network[!network$node1 %in% covariates, ]
    network <- network[!network$node2 %in% covariates, ]
    metadata <- get_metadata_for_columns(object=object, which_data=which_data, columns=cols_for_meta, 
                 names=c("name", "pathway"), index_of_names=rep("id", each=length(which_data)))
    out <- get_make_plotter_object(data=ggm_data, metadata=metadata, calc_type="genenet_ggm", 
      calc_info=paste("Cross-sectional GeneNet GGM results for:", paste(which_data, collapse=" & "), "at", timepoint, sep=" "), plot_type="network",
      style="visNetwork")
    return(out)

  })
#lol <- calc_ggm_genenet_crosssectional(object=data, which_data="lipid_data", threshold="li", all=FALSE, cols_for_meta=list(c("id", "sub_pathway")), timepoint="t0", covariates=c("Total_C", "HDL_C", "Total_TG", "VLDL_TG", "LDL_TG", "HDL_TG", "Age", "PTGENDER", "BMI"))


