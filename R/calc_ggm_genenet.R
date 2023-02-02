
#' An automated fucntion to calculate GGM from genenet crosssectional version
#' @description automated funtion that can be applied on metime_analyser object to obtain geneNet network along with threshold used
#' @param object S4 object of cĺass metime_analyser
#' @param which_data a character or a character vector naming the datasets of interest
#' @param threshold type of threshold to be used for extracting significant edges. 
#'      allowed inputs are "li", "FDR", "bonferroni"
#' @param all Logical to extract all edges without any pval correction
#' @param cols_for_meta list of character vector for extracting metadata of metabolites for plotting
#' @param covariates covariates to be used for this analysis
#' @param stratifications List to stratify data into a subset. Usage list(name=value)
#' @param name Name of the result
#' @param ... additional arguments for GeneNet
#' @return Network data as a plotter object
#' @export
setGeneric("calc_ggm_genenet", function(object, which_data, threshold, all, cols_for_meta, covariates, stratifications, name, ...) standardGeneric("calc_ggm_genenet"))
setMethod("calc_ggm_genenet", "metime_analyser", function(object, which_data, threshold, all, cols_for_meta, covariates, stratifications, name, ...) {
  
        if(is.null(stratifications)) {
          times <- object@list_of_row_data[[which_data[1]]]$time %>% unique()
        } else {
          times <- stratifications$time
        }
        data_list <- get_stratified_data(object=object, which_data=which_data,
                    stratifications=stratifications)
        data <- data_list[["data"]]
        row_data <- data_list[["row_data"]]

        if(length(times)>1) {
          data <- na.omit(data)
          data$subject <- rownames(data) %>% gsub(pattern="_[a-z|A-Z][0-9]+", replacement="")
          data$time <- rownames(data) %>% gsub(pattern="[a-z|A-Z][0-9]+_", replacement="")
          my_subjects <- data %>%     
              dplyr::select(all_of(c("subject","time"))) %>% 
              dplyr::filter(get("time") %in% stratifications$time)%>% 
              dplyr::group_by(subject) %>% 
              dplyr::count(subject) %>% 
              dplyr::filter(n == length(stratifications$time))
                  
          data <- data %>% 
            dplyr::filter(time %in% timepoints,
            subject %in% my_subjects[["subject"]])
          rm_col = intersect(names(data), c("adni_id","RID","rid","time","tp","subject", "id"))
          vars <- data %>% select(-c(all_of(rm_col))) %>% names()     
          # get full data
          data <- data %>% dplyr::arrange(time, subject) 
          n_subject = unique(data$subject) %>% length() %>% as.numeric()
          name_tp = as.numeric(unlist(lapply(strsplit(unique(data$time), split="t"), function(x) return(x[2]))))
          data <- longitudinal::as.longitudinal(x=as.matrix(data[,vars]), repeats=n_subject, time=name_tp)
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
            } else if(threshold=="li") {
              data <- data %>% as.matrix() %>% .[,] %>% as.data.frame()  
              cordat <- cor(data)
              eigenvals <- eigen(cordat)$values
              li.thresh <- sum( as.numeric(eigenvals >= 1) + (eigenvals - floor(eigenvals)) )
              met.ggm.edges.filtered <- met.ggm.edges[which(met.ggm.edges$pval < 0.05/li.thresh),]
            } else if(threshold=="bonferroni") {
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
          network <- met.ggm.edges.filtered[-edge2rem,]
          network <- network[!network$node1 %in% covariates, ]
          network <- network[!network$node2 %in% covariates, ]
          metadata <- get_metadata_for_columns(object=object, which_data=which_data, columns=cols_for_meta, 
                 names=c("name", "pathway"), index_of_names="id")
         out <- get_make_results(object=object, data=list(network), metadata=metadata, calc_type="genenet_ggm", 
              calc_info=paste("GeneNet GGM results for:", paste(which_data, collapse=" & "), "with",  
              ifelse(length(stratifications)>=1, paste(stratifications, collapse="_"), "full data"), sep=" "),
              name=name)
          out <- add_function_info(object=out, function_name="calc_ggm_genenet", 
            params=list(which_data=which_data, threshold=threshold, all=all,
              cols_for_meta=cols_for_meta, covariates=covariates, stratifications=stratifications, 
              name=name, ...))
          return(out)

        } else {
          
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
          out <- get_make_results(object=object, data=list(ggm_data), metadata=metadata, calc_type="genenet_ggm", 
              calc_info=paste("GeneNet GGM results for:", paste(which_data, collapse=" & "), "with",  
              ifelse(length(stratifications)>=1, paste(stratifications, collapse="_"), "full data"), sep=" "),
              name=name)
          out <- add_function_info(object=out, function_name="calc_ggm_genenet", 
            params=list(which_data=which_data, threshold=threshold, all=all,
              cols_for_meta=cols_for_meta, covariates=covariates, stratifications=stratifications, 
              name=name, ...))
          return(out)
        }

        

  })
#lol <- calc_ggm_genenet_crosssectional(object=data, which_data="lipid_data", threshold="li", all=FALSE, cols_for_meta=list(c("id", "sub_pathway")), timepoint="t0", covariates=c("Total_C", "HDL_C", "Total_TG", "VLDL_TG", "LDL_TG", "HDL_TG", "Age", "PTGENDER", "BMI"))


