#creating reference metime-analyser class that creates an object with full data

#' Constructor to generate an object of class metime_analyser. 
#' contains slots - list_of_data: For the list of all data matrices.
#'          - list_of_col_data: list of all the col data files in the same order.
#'          - list_of_row_data: list of all the row data files in the same order.
#'          - annotations: list with phenotype and medication. Each of which is character that represents 
#'                  the name of the aforementioned dataset types.   
#'  
#' @rdname metime_analyser
#' @export 
setClass("metime_analyser", slots=list(list_of_data="list", list_of_col_data="list", list_of_row_data="list", 
                 annotations="list")) 


#' Function to calculate dependent variables
#' @description An S4 method to be applied on the metime_analyser object so as to calculate dependent variables
#' @param object An object of class metime_analyser
#' @param which_x Name of the dataset to be used for training
#' @param which_y Name of the dataset to be used for testing
#' @param output_loc path to the parent directory where in the out file wíll be stored
#' @param file_name name of the out file
#' @param verbose Information provided on steps being processed
#' @return List of conservation index results
#' @export
#' 
#' 
setGeneric("calc_featureselection_boruta", function(object, which_x, which_y, verbose, output_loc, file_name) standardGeneric("calc_featureselection_boruta"))
setMethod("calc_featureselection_boruta", "metime_analyser", function(object, which_x,which_y, verbose=F, output_loc=getwd(), file_name="boruta") {
  
  # validate arguments
  #stopifnot(length(which_x)==1, length(which_y)==1, all(c(which_x,which_y) %in% object@list_of_data))
  
  x_data <- object@list_of_data[[which_x]]
  y_data <- object@list_of_data[[which_y]]
  
  x_data <- x_data[intersect(rownames(x_data), rownames(y_data)),]
  y_data <- y_data[rownames(x_data),]
  
  
  for(i in sample(colnames(y_data))){
    if(verbose) cat(i, "; ")
    if(paste0(file_name, ".rds") %in% list.files(output_loc)) my_results=readRDS(file=paste0(output_loc, "/", file_name, ".rds"))
      else my_results = data.frame()
      
    my_model = Boruta::Boruta(x=x_data, y=y_data[,i],
                              pValue=0.01,
                              doTrace=0,
                              mcAdj=TRUE)
    my_stats=Boruta::attStats(my_model) %>%
      as.data.frame() %>% 
      dplyr::mutate(id_x=i,
                    id_y=rownames(.[]))
    
    saveRDS(object=rbind(my_results, my_stats),file=paste0(output_loc, "/", file_name, ".rds")) 
    
  }
})



#' Function to calculate metabotype conservation index
#' @description Method applied on the object metime_analyser to calculate the metabotype conservation index
#' @examples #calculating metabotype_conservation_index 
#' out <- calc_metabotype_conservation(object=metime_analyser_object, which_data="Name of the dataset")
#' @param object An object of class metime_analyser
#' @param which_data Name of the dataset to be used
#' @param timepoints character vector with timepoints of interest
#' @param verbose Information provided on steps being processed
#' @return List of conservation index results
#' @export

setGeneric("calc_conservation_metabotype", function(object, which_data, timepoints, verbose) standardGeneric("calc_conservation_metabotype"))
setMethod("calc_conservation_metabotype", "metime_analyser", function(object, which_data, timepoints, verbose=F) {
  #define data to be processed
  #data_position <- which(names(object@list_of_data) %in% which_data)
  out=list()
  for(i in which_data){
    data_position <- which(names(object@list_of_data) %in% i)
    data_var <- object@list_of_col_data[[data_position]]$id
    
    data_merged <- object@list_of_data[[data_position]] %>% 
      dplyr::mutate(id = rownames(.[])) %>% 
      dplyr::left_join(
        object@list_of_row_data[[data_position]] %>% 
          dplyr::mutate(timepoint = as.numeric(timepoint)) %>% 
          dplyr::select(id, timepoint, rid),
        by="id"
      )
    
    tp_split <- utils::combn(x = unique(data_merged$timepoint), m=2, simplify = T) %>% 
      base::t() %>% 
      base::as.data.frame()
    timepoints <- unlist(lapply(strsplit(timepoints, split="t"), function(x) return(x[2])))
    tp_split <- tp_split[which(as.character(tp_split$V1) %in% timepoints & as.character(tp_split$V2) %in% timepoints), ]
    
    id_split <- data_merged %>% select(id, rid, timepoint) %>% 
      tidyr::spread(timepoint, id)
    
    
    if(verbose) cat("Processing dataframe", i,": ")

    out[[i]] = lapply(1:nrow(tp_split), function(x) {
      
      if(verbose) cat(tp_split[x,] %>% paste0(collapse="vs"), "; ")

      v1_data <- data_merged %>% 
                dplyr::filter(id %in% id_split[,paste0(tp_split$V1[[x]])]) %>% 
                 dplyr::select(-id, -timepoint,-rid) %>% t() %>% as.data.frame()

      colnames(v1_data) = intersect(id_split[,paste0(tp_split$V1[[x]])], data_merged$id)
      
      v2_data <- data_merged %>% dplyr::filter(id %in% id_split[,paste0(tp_split$V2[[x]])]) %>% 
        dplyr::select(-id, -timepoint,-rid) %>% t() %>% as.data.frame()
      colnames(v2_data) = intersect(id_split[,paste0(tp_split$V2[[x]])], data_merged$id)
      
      cor_mat = cbind(v1_data,v2_data) %>% 
        cor(use="pairwise.complete.obs")
      
      ci_out <- lapply(1:nrow(id_split), function(y) {
        v2_id = id_split[,paste0(tp_split$V2[[x]])][y]
        v1_id = id_split[,paste0(tp_split$V1[[x]])][y]
        if(all(!is.na(c(v2_id, v1_id)))) {
          
          ci_info <- cor_mat[intersect(id_split[,paste0(tp_split$V2[[x]])], data_merged$id),v1_id] %>% sort(decreasing = T)
          rank <- which(names(ci_info)==v2_id)
          return(data.frame(
            id=v1_id,
            from_tp=tp_split[x,"V1"],
            to_tp=tp_split[x,"V2"],
            var=id_split$rid[which(id_split[[paste0(tp_split$V2[[x]])]]==v2_id)],
            nsubject=length(ci_info),
            rank=ifelse(length(rank)>0, rank, NA),
            stringsAsFactors = F))
        } 
      }) %>% 
        do.call(what=rbind.data.frame)
      
      ci_out = ci_out %>%
        dplyr::mutate(ci=(1 - ((rank -1)/(nsubject-1)))) %>% 
        dplyr::arrange(ci) %>% 
        dplyr::mutate(y=ci, 
                      x=1:nrow(.[])) %>% 
        dplyr::select(x,y,ci,id, from_tp, to_tp, nsubject,rank)
      rownames(ci_out) <- ci_out$id
      return(ci_out)
    })
    
    names(out[[i]]) = paste0(tp_split$V1,"vs", tp_split$V2)
  }
  
  return(out)
})


#' Function to calculate metabolite conservation index
#' @description Method applied on the object metime_analyser to calculate the metabotype conservation index
#' @examples #calculating metabolite_conservation_index 
#' out <- calc_metabolite_conservation(object=metime_analyser_object, which_data="Name of the dataset")
#' @param object An object of class metime_analyser
#' @param which_data Name of the dataset to be used
#' @param timepoints character vector with timepoints of interest
#' @param verbose Information provided on steps being processed
#' @return List of conservation index results
#' @export
#' 
#' 
setGeneric("calc_conservation_metabolite", function(object, which_data, timepoints, verbose) standardGeneric("calc_conservation_metabolite"))
setMethod("calc_conservation_metabolite", "metime_analyser", function(object, which_data, timepoints, verbose=F) {
  #define data to be processed
  #data_position <- which(names(object@list_of_data) %in% which_data)
  out=list()
  for(i in which_data) {
    data_position <- which(names(object@list_of_data) %in% i)
    data_var <- object@list_of_col_data[[data_position]]$id
    
    data_merged <- object@list_of_data[[data_position]] %>% 
      dplyr::mutate(id = rownames(.[])) %>% 
      dplyr::left_join(object@list_of_row_data[[data_position]] %>% 
          dplyr::mutate(timepoint = as.numeric(timepoint)) %>% 
          dplyr::select(id, timepoint, rid), 
          by="id")
    
    tp_split <- utils::combn(x = unique(data_merged$timepoint), m=2, simplify = T) %>% 
      base::t() %>% 
      base::as.data.frame()

    timepoints <- unlist(lapply(strsplit(timepoints, split="t"), function(x) return(x[2])))
    tp_split <- tp_split[which(as.character(tp_split$V1) %in% timepoints & as.character(tp_split$V2) %in% timepoints), ]
    
    id_split <- data.frame(met=rep(data_var, each=length(unique(data_merged$timepoint))),
      timepoint=rep(unique(data_merged$timepoint), length(data_var)),
      stringsAsFactors = F) %>% 
        dplyr::mutate(id=paste0(met, "_",timepoint)) %>% 
          tidyr::spread(timepoint, id)
    
    
    if(verbose) cat("Processing dataframe", i,": ")
      out[[i]] = lapply(1:nrow(tp_split), function(x){
        
        if(verbose) cat(tp_split[x,] %>% paste0(collapse="vs"), "; ")
        
        

        v1_data <- data_merged %>% dplyr::filter(timepoint == tp_split[x,"V1"]) %>% 
                    dplyr::select(-id, -timepoint) 
        colnames(v1_data)[-length(colnames(v1_data))]  = paste0(names(v1_data)[-length(names(v1_data))], "_", tp_split$V1[x])
        
        v2_data <- data_merged %>% dplyr::filter(timepoint == tp_split[x,"V2"])  %>% 
                    dplyr::select(-id, -timepoint)
        colnames(v2_data)[-length(colnames(v2_data))]  = paste0(names(v2_data)[-length(names(v2_data))], "_",tp_split$V2[x])

        if(length(v1_data$rid) > length(v2_data$rid)) {
             v1_data <- v1_data[v1_data$rid %in% v2_data$rid, ]
             v2_data <- v2_data[v2_data$rid %in% v1_data$rid, ]
        } else {
              v2_data <- v2_data[v2_data$rid %in% v1_data$rid, ]
              v1_data <- v1_data[v1_data$rid %in% v2_data$rid, ] 
        }

        v1_data <- v1_data %>% dplyr::select(-rid)
        v2_data <- v2_data %>% dplyr::select(-rid)

        cor_mat = cbind(v1_data,v2_data) %>% 
          cor(use="pairwise.complete.obs")
        
        ci_out <- lapply(1:nrow(id_split), function(y){
          v2_id = id_split[,paste0(tp_split$V2[[x]])][y]
          v1_id = id_split[,paste0(tp_split$V1[[x]])][y]
          if(all(!is.na(c(v2_id, v1_id)))){
            ci_info <- cor_mat[intersect(rownames(cor_mat),id_split[,paste0(tp_split$V2[[x]])]),v1_id] %>% sort(decreasing = T)
            rank <- which(names(ci_info)==v2_id)
            
            data.frame(
              id=id_split$met[y],
              from_tp=tp_split[x,"V1"],
              to_tp=tp_split[x,"V2"],
              nsubject=length(ci_info),
              rank=ifelse(length(rank)>0, rank, NA),
              stringsAsFactors = F)
          }
          
        }) %>% 
          do.call(what=rbind.data.frame)
        
        ci_out = ci_out %>%
          dplyr::mutate(ci=(1 - ((rank -1)/(nsubject-1)))) %>% 
          dplyr::arrange(ci) %>% 
          dplyr::mutate(y=ci, 
                        x=1:nrow(.[])) %>% 
          dplyr::select(x,y,ci,id, from_tp, to_tp, nsubject,rank)
        rownames(ci_out) <- ci_out$id
        return(ci_out)
      })
      
      names(out[[i]]) = paste0(tp_split$V1,"vs", tp_split$V2)
    }
  
  return(out)
})


#' Function to calculate dimensionality reduction methods such as tsne, umap and pca.
#' @description A method to apply on s4 object of class metime_analyse in order to obtain information after dimensionality reduction on a dataset/s
#' @examples
#' #calculate PCA
#' pca <- calc_dimensionality_reduction(object=metime_analyser_object, which_data="name/s of the dataset/s", type="PCA")
#' #calculate UMAP
#' pca <- calc_dimensionality_reduction(object=metime_analyser_object, which_data="name/s of the dataset/s", type="UMAP")
#' #calculate tSNE
#' pca <- calc_dimensionality_reduction(object=metime_analyser_object, which_data="name/s of the dataset/s", type="tSNE")
#' @param object An object of class metime_analyser
#' @param which_data a character vector - Names of the dataset from which the samples will be extracted
#' @param type type of the dimensionality reduction method to be applied. Accepted inputs are "UMAP", "tSNE", "PCA"
#' @param ... additional arguments that can be passed on to prcomp(), M3C::tsne() and umap::umap()
#' @return  a list with two dataframes containing the dimensionality reduction information 1) samples - data of the individuals(".$samples")
#'                     2) metabs - data of the metabolites(".$metabs")
#' @export
setGeneric("calc_dimensionality_reduction", function(object, which_data, type, ...) standardGeneric("calc_dimensionality_reduction"))

setMethod("calc_dimensionality_reduction", "metime_analyser", function(object, which_data, type, ...) {
      if(length(which_data) > 1) {
        object <- mod_extract_common_samples(object)
        data <- object@list_of_data[names(object@list_of_data) %in% which_data]
        data <- lapply(data, function(x) {
            return(x[sort(rownames(x)),])
          })
        data <- as.data.frame(do.call(cbind, data))
        colnames(data) <- unlist(lapply(strsplit(colnames(data), split="_data.", fixed=TRUE), function(x) return(x[2])))
      } else {
        data <- object@list_of_data[names(object@list_of_data) %in% which_data][[1]]
      }
      rm_col = intersect(names(data), c("adni_id","RID","rid","timepoint","tp","subject", "id"))
      data <- data %>% select(-c(rm_col))
      data <- na.omit(data)
      if(type %in% "PCA") {
        pca_metabs <- prcomp(t(data), scale.=T, center=T, ...)
        pca_individuals <- prcomp(data, scale.=T, center=T, ...)
        dr_data_metabs <- as.data.frame(pca_metabs$x[,1:2])
        dr_data_samples <- as.data.frame(pca_individuals$x[,1:2])
      } else if(type %in% "UMAP") {
        umap_individuals <- umap::umap(data, ...)
        dr_data_samples <- as.data.frame(umap_individuals$layout)
        colnames(dr_data_samples) <- c("UMAP1", "UMAP2")
        umap_metabs <- umap::umap(t(data), ...)
        dr_data_metabs <- as.data.frame(umap_metabs$layout)
        colnames(dr_data_metabs) <- c("UMAP1", "UMAP2")
      } else if(type %in% "tSNE") {
        tsne_samples <- tsne(t(data), ...)
        tsne_metabs <- tsne(data, ...)
        dr_data_metabs <- as.data.frame(tsne_metabs$data)
        rownames(dr_data_metabs) <- colnames(data)
        dr_data_samples <- as.data.frame(tsne_samples$data)
        rownames(dr_data_samples) <- rownames(data)
      }
      return(list(metabs=dr_data_metabs, samples=dr_data_samples))
  })


#' Function to calculate correlation 
#'
#' @description calculate pairwise correlations
#' This function creates a dataframe for plotting from a dataset.
#' @examples # Example to calculate correlations
#' dist <- calc_correlation(object=metime_analyser_object, which_data="name of the dataset", 
#'           method="pearson")
#' @param object S4 Object of class metime_analyser
#' @param which_data specify datasets to calculate on. One or more possible
#' @param method default setting: method="pearson", Alternative "spearman" also possible
#' @return data.frame with pairwise results
#' @export
setGeneric("calc_correlation_pairwise", function(object, which_data, method) standardGeneric("calc_correlation_pairwise"))
setMethod("calc_correlation_pairwise", "metime_analyser", function(object, which_data, method="pearson"){
  stopifnot(all(which_data %in% names(object@list_of_data)))
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    return(data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    ))
  }
    
  my_data <-  lapply(which_data, function(x) object@list_of_data[[x]] %>% 
                         dplyr::mutate(id=rownames(.[]))) %>% 
      plyr::join_all(by="id", type="inner") %>% 
      `rownames<-`(.[,"id"]) %>% 
      dplyr::select(-id)
    
    # calculate correlation matrix and pvalues
    cor_mat <- my_data %>% 
      as.matrix() %>% 
      Hmisc::rcorr(type=method)
    out=flattenCorrMatrix(cor_mat$r, cor_mat$P) %>% 
      dplyr::mutate(type="cor") %>% 
      dplyr::rename("dist"="cor", "cut_p"="p") 
    return(out)
})

#' Function to calculate dissimilarity using distance measures 
#'
#' @description calculate pairwise distances
#' This function creates a dataframe for plotting from a dataset.
#' @examples # Example to calculate pairwise distances
#' dist <- calc_pairwise_distance(object=metime_analyser_object, which_data="name of the dataset", 
#'           method="euclidean")
#' @param object S4 Object of class metime_analyser
#' @param which_data specify datasets to calculate on. One or more possible
#' @param method default setting: method="euclidean", Alternative "maximum","minimum",
#' "manhattan","canberra","minkowski" are also possible
#' @return data.frame with pairwise results
#' @export

setGeneric("calc_distance_pairwise", function(object, which_data, method) standardGeneric("calc_distance_pairwise"))
setMethod("calc_distance_pairwise", "metime_analyser", function(object, which_data, method="euclidean") {
  stopifnot(all(which_data %in% names(object@list_of_data)))
  
  flattenCorrMatrix <- function(cormat) {
    ut <- upper.tri(cormat)
    return(data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      dist  =(cormat)[ut]
    ))
  }
  
  my_data <-  lapply(which_data, function(x) object@list_of_data[[x]] %>% 
                       dplyr::mutate(id=rownames(.[]))) %>% 
    plyr::join_all(by="id", type="inner") %>% 
    `rownames<-`(.[,"id"]) %>% 
    dplyr::select(-id)
  
  if(dist %in% c("euclidean","maximum","minimum","manhattan","canberra","minkowski")){
    
    out <- my_data %>%
      stats::dist(method = method) %>%
      as.matrix() %>% 
      as.data.frame() %>% 
      flattenCorrMatrix() %>% 
      dplyr::mutate(type=method)
  }
  else{
    out=NA
  }
  return(out)
})


#' Function to calculate students t-test
#' @description Method for S4 object of class metime_analyser for performing t-test
#' @param object S4 object of class metime_analyser
#' @param which_data dataset or datasets to be used for the analysis
#' @param timepoints two timepoints of interest to perform the test on
#' @param split_var split variable for testing such as diagnostic group etc
#' @return t-test result as a list or a list of t-test results 
#' @export
setGeneric("calc_ttest", function(object, which_data, timepoints, split_var) standardGeneric("calc_ttest"))
setMethod("calc_ttest", "metime_analyser", function(object, which_data, timepoints, split_var) {
        if(length(which_data) > 1) object <- mod_extract_common_samples(object)
        list_of_data <- mod_split_acc_to_time(object)
        list_of_data <- list_of_data[names(list_of_data) %in% which_data]
        list_of_data <- lapply(list_of_data, function(x) {
              x <- x[names(x) %in% timepoints]
              count <- 1
              x <- lapply(x, function(a) {
                    a <- a[ ,order(colnames(a))]
                    colnames(a) <- paste(colnames(a), "_time:", timepoints[count], sep="")
                    count <- count + 1
                    return(a)
              })
              return(x)
        })
        if(!is.null(split_var)) {
            for(i in 1:length(list_of_data)) {

            }
        }
        

  })


#' An automated fucntion to calculate GGM from genenet longitudnal version
#' @description automated funtion that can be applied on metime_analyser object to obtain geneNet network along with threshold used
#' @param object S4 object of cĺass metime_analyser
#' @param which_data a character or a character vector naming the datasets of interest
#' @param threshold type of threshold to be used for extracting significant edges
#' @param timepoints timepoints of interest that are to be used to build networks(as per timepoints in rows)
#' @param all Logical to get all edges without any cutoff.
#' @param ... addtional arguments for genenet network
#' @return Network data with edgelist, partial correlation values and associated p-values and corrected p-values 
#' @export
setGeneric("calc_ggm_genenet_longitudnal", function(object, which_data, threshold, timepoints, all, ...) standardGeneric("calc_ggm_genenet_longitudnal"))
setMethod("calc_ggm_genenet_longitudnal", "metime_analyser", function(object, which_data, threshold, timepoints, all, ...) {
    #sanity checks
    stopifnot(timepoints %in% c("t0","t12","t24"))
    stopifnot(threshold %in% c("li","bonferroni","FDR", NULL))
    #Extracting data that is needed
    object@list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
    object@list_of_row_data <- object@list_of_row_data[names(object@list_of_row_data) %in% which_data]
    #converting prepped data to full ggm network
    if(length(which_data > 1)) {
        object <- mod_extract_common_samples(object)
        list_of_data <- unname(object@list_of_data)
        data <- do.call(cbind, list_of_data)
    } else {
        data <- object@list_of_data[[1]]
    } 
    data <- na.omit(data)
    data$subject <- unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[1])))
    data$timepoint <- unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[2])))
    my_rid <- data %>%     
              dplyr::select(all_of(c("subject","timepoint"))) %>% 
              dplyr::filter(get("timepoint") %in% timepoints)%>% 
              dplyr::group_by(subject) %>% 
              dplyr::count(subject) %>% 
              dplyr::filter(n == length(timepoints))
                  
    data <- data %>% 
            dplyr::filter(timepoint %in% timepoints,
            subject %in% my_rid[["subject"]])
    rm_col = intersect(names(data), c("adni_id","RID","rid","timepoint","tp","subject", "id"))
    vars <- data %>% select(-c(all_of(rm_col))) %>% names()


    # get full data
    data <- data %>% dplyr::arrange(timepoint, subject) 
    n_subject = unique(data$subject) %>% length() %>% as.numeric()
    name_tp = as.numeric(unlist(lapply(strsplit(unique(data$timepoint), split="t"), function(x) return(x[2]))))

    data <- longitudinal::as.longitudinal(x=as.matrix(data[,vars]), repeats=n_subject, time=name_tp)

    network <- get_ggm_genenet(data=data, threshold=threshold, all=all, ...)
       
    return(network)
})


#' An automated fucntion to calculate GGM from multibipartite lasso approach
#' @description automated funtion that can be applied on s4 object of class metime_analyser to calculate a network using
#' multibipartite lasso
#' @param object S4 object of cĺass metime_analyser
#' @param which_data a character or a character vector naming the datasets of interest
#' @param alpha tuning parameter for lasso + ridge regression in glmnet
#' @param nfolds nfolds for cv.glmnet
#' @param timepoints timepoints of interest that are to be used to build networks(as per timepoints in rows)
#' @return Network data with edges and their respective betas
#' @export
setGeneric("calc_ggm_multibipartite_lasso", function(object, which_data, alpha, nfolds, timepoints) standardGeneric("calc_ggm_multibipartite_lasso"))
setMethod("calc_ggm_multibipartite_lasso", "metime_analyser", function(object, which_data, alpha, nfolds, timepoints) {
          if(length(which_data) > 1) object <- mod_extract_common_samples(object)
          list_of_data <- mod_split_acc_to_time(object)
          list_of_data <- list_of_data[names(list_of_data) %in% which_data]
          list_of_data <- lapply(list_of_data, function(x) return(x[names(x) %in% timepoints]))
          count <- 1
          edge_lists <- list()
          for(i in 1:length(timepoints)) {
              list_of_mats <- list()
              list_of_mats <- lapply(list_of_data, function(x) {
                            x <- x[names(x) %in% timepoints[i]]
                            x <- lapply(x, function(a) return(a[order(rownames(a)), ]))
                            x <- lapply(x, function(a) return(as.matrix(a)))
                            return(x[[1]]) 
              })
              names(list_of_mats) <- which_data
              #glmnet results are stored here in a list
              list <- get_betas_for_multibipartite_lasso(list_of_mats=list_of_mats, 
                                            alpha=alpha, nfolds=nfolds)
              #list to check and collect all the edges that have non-zero coefficient 
              list_to_check <- list()
              #list here is the list with all the details from glmnet
              for(j in 1:length(list)) {
                  #creating a list to store the edges
                  list_of_edges <- list()
                  #list[[1]] is the data for first metabolite
                  for(k in 1:length(list[[j]])) {
                      #coeffs is the vector with all the coefficient values
                      coeffs <- coef(list[[j]][[k]])[,1]
                      #removing zero coefficients
                      coeffs <- coeffs[!(coeffs==0)]
                      #removing the intercept data
                      coeffs <- coeffs[-1]
                      #source is the metab names
                      source <- rep(names(list[[j]])[k], each=length(coeffs))
                      #target is the name of the metabolites that have non-zero coeffs
                      target <- names(coeffs)
                      #storing data as a nested list
                      list_of_edges[[k]] <- cbind(source, target, coeffs) 
                      names(list_of_edges)[k] <- source[1]
                     
                }
                #list_to_check will now store the edges
                list_to_check[[j]] <- list_of_edges
              }
            edge_list <- unlist(list_to_check, recursive=FALSE)
            edge_list <- as.data.frame(do.call(rbind, edge_list)) 
            uids <- c()
            #ERROR FOUND HERE AND CORRECTED!!!!!!!!
            for(m in 1:length(edge_list[,1])) {
                vec <- c()
                vec[1] <- as.character(edge_list$source[m])
                vec[2] <- as.character(edge_list$target[m])
                vec <- vec[order(vec)]
                uids[m] <- paste(vec, collapse="_")
            }
            edge_list <- cbind(edge_list, uids)
            edge_list <- as.data.frame(edge_list)
            colnames(edge_list) <- c("node1", "node2", "coeffs", "uids")
            check <- edge_list[,c("uids", "coeffs")]
            check$uids <- as.character(check$uids)
            check$coeffs <- as.numeric(as.character(check$coeffs))
            check <- reshape(transform(check, time=ave(coeffs, uids, FUN=seq_along)), idvar="uids", direction="wide")
            final_edge_list <- na.omit(check)
            final_edge_list$uids <- as.character(final_edge_list$uids)
            final_edge_list$node1 <- unlist(lapply(strsplit(final_edge_list$uids, split="_"), function(x) return(x[1])))
            final_edge_list$node2 <- unlist(lapply(strsplit(final_edge_list$uids, split="_"), function(x) return(x[2])))
            final_edge_list$uids <- NULL
            edge_lists[[i]] <- final_edge_list            
          }
          names(edge_lists) <- timepoints
          return(edge_lists)
  })


#lol <- calc_ggm_multibipartite_lasso(object=data, alpha=1, which_data=c("lipid_data", "nmr_data"), nfolds=3, timepoints="t0")

#' An automated function to caluclate temporal network with lagged model
#' @description calculates temporal networks for each dataset with a lagged model as used in graphical VAR
#' @param object S4 object of class metab_analyser
#' @param lag which lagged model to use. 1 means one-lagged model, similary 2,3,..etc
#' @param which_data dataset or datasets to be used
#' @param timepoints timepoints of interest that are to be used to build networks(in the order of measurement)
#' @param alpha parameter for regression coefficient
#' @param nfolds nfolds parameter for glmnet style of regression
#' @return temporal network data with edgelist and regression values
#' @export
setGeneric("calc_temporal_ggm", function(object, which_data, lag, timepoints, alpha, nfolds) standardGeneric("calc_temporal_ggm"))
setMethod("calc_temporal_ggm", "metime_analyser", function(object, which_data, lag, timepoints, alpha, nfolds) {
        if(length(which_data) > 1) object <- mod_extract_common_samples(object)
        object@list_of_data <- mod_split_acc_to_time(object)
        list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
        list_of_data <- lapply(list_of_data, function(x) {
              x <- x[names(x) %in% timepoints]
              x <- lapply(x, function(a) return(a[order(rownames(a)), ]))
              return(x)
          })
        final_data <- list()
        list_of_data <- unname(list_of_data)
        for(i in 1:length(timepoints)) {
            timepoint_list <- lapply(list_of_data, function(x) return(x[[i]]))
            final_data[[i]] <- do.call(cbind, timepoint_list)
        }
        names(final_data) <- timepoints
        count <- 1
        model_seqs <- list()
        while(lag + count <= length(timepoints)) {
            model_seqs[[count]] <- timepoints[seq(from=lag+count, to=count, by=-1)]
            count <- count + 1
        }
        models <- list()
        for(i in 1:length(model_seqs)) {
            network_data <- final_data[names(final_data) %in% model_seqs[[i]]]
            count <- 1
            network_data <- lapply(network_data, function(x) {
                    colnames(x) <- paste(colnames(x), model_seqs[[i]][count], sep="_time:")
                    rownames(x) <- unlist(lapply(strsplit(rownames(x), split="_"), function(x) return(x[1])))
                    count <<- count +1 
                    return(x) 
              })
            network_data <- unname(network_data)
            list_of_names <- lapply(network_data, function(x) {
                return(rownames(x))
            })
            common_samples <- Reduce(intersect, list_of_names)
            network_data <- lapply(network_data, function(x) {
                x <- x[rownames(x) %in% common_samples, ]
                x <- x[order(rownames(x)), ]
                return(x)
            })
            network_data <- do.call(cbind, network_data)
            ymat <- network_data[ ,grep(model_seqs[[i]][1] ,colnames(network_data), value=TRUE)]
            xmat <- network_data[ ,!(colnames(network_data) %in% colnames(ymat))]
            ymat <- as.matrix(na.omit(ymat))
            xmat <- as.matrix(na.omit(xmat))
            fit_list <- list()
            for(k in 1:ncol(ymat)) {
                fit_list[[k]] <- cv.glmnet(x=xmat, y=ymat[,k], alpha=alpha, nfolds=nfolds)
                coeffs <- coef(fit_list[[k]])[,1]
                coeffs <- coeffs[!(coeffs==0)]
                coeffs <- coeffs[-1]
                target <- names(coeffs)
                source <- rep(colnames(ymat)[k], each=length(coeffs))
                fit_list[[k]] <- as.data.frame(cbind(source, target, coeffs))
            }
            fit_list <- unname(fit_list)
            models[[i]] <- as.data.frame(do.call(rbind, fit_list))
            names(models)[i] <- paste(model_seqs[[i]], collapse="-")
            colnames(models[[i]]) <- c("node1", "node2", "coeffs")
            models[[i]]$label <- paste(unlist(lapply(strsplit(as.character(models[[i]]$node1), split="_time:"), function(x) return(x[2]))),
                                      unlist(lapply(strsplit(as.character(models[[i]]$node2), split="_time:"), function(x) return(x[2]))),
                                        sep="-")
            models[[i]]$node1 <- unlist(lapply(strsplit(as.character(models[[i]]$node1), split="_time:"), function(x) return(x[1])))
            models[[i]]$node2 <- unlist(lapply(strsplit(as.character(models[[i]]$node2), split="_time:"), function(x) return(x[1])))
            
        }
        return(models)
  })


#' An automated fucntion to calculate GGM from genenet crosssectional version
#' @description automated funtion that can be applied on metime_analyser object to obtain geneNet network along with threshold used
#' @param object S4 object of cĺass metime_analyser
#' @param which_data a character or a character vector naming the datasets of interest
#' @param threshold type of threshold to be used for extracting significant edges. 
#'      allowed inputs are "li", "FDR", "bonferroni"
#' @param timepoints timepoints of interest that are to be used to build networks(as per timepoints in rows)
#' @param all Logical to extract all edges without any pval correction
#' @param ... additional arguments for GeneNet
#' @return Network data with edgelist, partial correlation values and associated p-values and corrected p-values 
#' @export
setGeneric("calc_ggm_genenet_crosssectional", function(object, which_data, threshold, timepoint, all, ...) standardGeneric("calc_ggm_genenet_crosssectional"))
setMethod("calc_ggm_genenet_crosssectional", "metime_analyser", function(object, which_data, threshold, timepoint, all, ...) {
        if(length(which_data) > 1) object <- mod_extract_common_samples(object)
        stopifnot(threshold %in% c("li","bonferroni","FDR", NULL))
        object@list_of_data <- mod_split_acc_to_time(object)
        list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
        list_of_data <- lapply(list_of_data, function(x) return(x[names(x) %in% timepoint]))
        list_of_data <- lapply(list_of_data, function(x) return(x[[1]][order(rownames(x[[1]])),]))
        list_of_data <- unname(list_of_data)
        data <- as.data.frame(do.call(cbind, list_of_data))
        rm_col = intersect(colnames(data), c("adni_id","RID","rid","timepoint","tp","subject", "id"))
        data <- data %>% select(-c(all_of(rm_col)))
        data <- na.omit(data)
        this_mat <- as.matrix(apply(data, 2, as.numeric))
        pcor_mat <- ggm.estimate.pcor(as.matrix(this_mat), method = "dynamic", verbose = F)
        # compute p-values of edges
        pval_mat <- network.test.edges(pcor_mat, plot = F, verbose = F)
        # p-value correction of edges
        pval_mat$p.adj.bh <- p.adjust(pval_mat$pval, method="BH")
        pval_mat$p.adj.bon <- p.adjust(pval_mat$pval, method="bonferroni")    
        ggm_thresh <- 0.05
        # extract edge list
        tmp <- pcor_mat %>% graph_from_adjacency_matrix(mode='undirected', weighted = T) %>% igraph::simplify()
        ggm_edges <- cbind.data.frame(get.edgelist(tmp), edge_attr(tmp)$weight)
        names(ggm_edges) <- c("node1", "node2", "pcor_val")
        if(all) {
           return(ggm_edges)
        } else {
          if(threshold %in% "bonferroni") {
            #filter edges based on p-values - bonferroni
            ggm_data <-  ggm_edges %>% filter(abs(pcor_val)>=min(abs(pval_mat$pcor[pval_mat$p.adj.bon<=ggm_thresh]))) 
          } else if(threshold %in% "FDR"){
           # filter edges based on p-values - BH 
            ggm_data <-  ggm_edges %>% filter(abs(pcor_val)>=min(abs(pval_mat$pcor[pval_mat$p.adj.bh<=ggm_thresh])))
          } else if(threshold %in% "li") {
            data <- this_mat %>% as.matrix() %>% .[,] %>% as.data.frame()  
            cordat <- cor(data)
            eigenvals <- eigen(cordat)$values
            li.thresh <- sum( as.numeric(eigenvals >= 1) + (eigenvals - floor(eigenvals)) )
            ggm_data <- ggm_edges %>% filter(abs(pcor_val)>=min(abs(pval_mat$pcor[pval_mat$pval<=0.05/li.thresh])))
          }
          return(ggm_data)        
        }

  })


##search for metabolic networks in ad
#Download any lipid based network from kegg
#write a fucntion to compare edges and check out elisas paper again
#see the outliers of pca - in tsne and umap
#remove unnecessary texts in rmarkdowns
#read parafac and repeated asca+ papers and see which one is suitable for us
#R packages that can be used to enrich networks for metabolites
#think of ideas for R shiny app, see humet and adatlas. test app to see if the backend needs to be changed