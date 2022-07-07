#creating reference metab-analyser class that creates an object with full data

#' Constructor to generate an object of class metab_analyser. 
#' contains slots - list_of_data: For the list of all data matrices.
#'          - list_of_col_data: list of all the col data files in the same order.
#'          - list_of_row_data: list of all the row data files in the same order.
#'          - annotations: list with phenotype and medication. Each of which is character that represents 
#'                  the name of the aforementioned dataset types.   
#'  
#' @rdname metab_analyser
#' @export 
setClass("metab_analyser", slots=list(list_of_data="list", list_of_col_data="list", list_of_row_data="list", 
                 annotations="list")) 


#' Function to calculate metabotype conservation index
#' @param object An object of class metab_analyser
#' @param which_data Name of the dataset to be used
#' @param verbose Information provided on steps being processed
#' @return List of conservation index results
#' @export
#' 
#' 
setGeneric("calc_boruta", function(object, which_x, which_y, verbose, output_loc, file_name) standardGeneric("calc_boruta"))
setMethod("calc_boruta", "metab_analyser", function(object, which_x,which_y, verbose=F, output_loc=getwd(), file_name="boruta") {
  
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
#' @param object An object of class metab_analyser
#' @param which_data Name of the dataset to be used
#' @param verbose Information provided on steps being processed
#' @return List of conservation index results
#' @export

setGeneric("calc_metabotype_conservation", function(object, which_data, verbose) standardGeneric("calc_metabotype_conservation"))
setMethod("calc_metabotype_conservation", "metab_analyser", function(object, which_data, verbose=F) {
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
    
    id_split <- data_merged %>% select(id, rid, timepoint) %>% 
      tidyr::spread(timepoint, id)
    
    
    if(verbose) cat("Processing dataframe", i,": ")
    out[[i]] = lapply(1:nrow(tp_split), function(x){
      
      if(verbose) cat(tp_split[x,] %>% paste0(collapse="vs"), "; ")


      v1_data <- data_merged %>% dplyr::filter(id %in% id_split[,paste0(tp_split$V1[[x]])]) %>% 
        dplyr::select(-id, -timepoint,-rid) %>% t() %>% as.data.frame()
      colnames(v1_data) = intersect(id_split[,paste0(tp_split$V1[[x]])], data_merged$id)
      
      v2_data <- data_merged %>% dplyr::filter(id %in% id_split[,paste0(tp_split$V2[[x]])]) %>% 
        dplyr::select(-id, -timepoint,-rid) %>% t() %>% as.data.frame()
      colnames(v2_data) = intersect(id_split[,paste0(tp_split$V2[[x]])], data_merged$id)
      
      cor_mat = cbind(v1_data,v2_data) %>% 
        cor(use="pairwise.complete.obs")
      
      ci_out <- lapply(1:nrow(id_split), function(y){
        v2_id = id_split[,paste0(tp_split$V2[[x]])][y]
        v1_id = id_split[,paste0(tp_split$V1[[x]])][y]
        if(all(!is.na(c(v2_id, v1_id)))){
          ci_info<-cor_mat[intersect(id_split[,paste0(tp_split$V2[[x]])], data_merged$id),v1_id] %>% sort(decreasing = T)
          rank<-which(names(ci_info)==v2_id)
          
          data.frame(
            id=v1_id,
            from_tp=tp_split[x,"V1"],
            to_tp=tp_split[x,"V2"],
            var=id_split$rid[which(id_split[[paste0(tp_split$V2[[x]])]]==v2_id)],
            nsubject=length(ci_info),
            rank=ifelse(length(rank)>0, rank, NA),
            stringsAsFactors = F
          )
        }else{

        }
          
        
      }) %>% 
        do.call(what=rbind.data.frame)
      
      ci_out = ci_out %>%
        dplyr::mutate(ci=(1 - ((rank -1)/(nsubject-1)))) %>% 
        dplyr::arrange(ci) %>% 
        dplyr::mutate(y=ci, 
                      x=1:nrow(.[])) %>% 
        dplyr::select(x,y,ci,id, from_tp, to_tp, nsubject,rank)
      
      ci_out
    })
    
    names(out[[i]]) = paste0(tp_split$V1,"vs", tp_split$V2)
  }
  
  return(out)
})


#' Function to calculate metabolite conservation index
#' @param object An object of class metab_analyser
#' @param which_data Name of the dataset to be used
#' @param verbose Information provided on steps being processed
#' @return List of conservation index results
#' @export
#' 
#' 
setGeneric("calc_metabolite_conservation", function(object, which_data, verbose) standardGeneric("calc_metabolite_conservation"))
setMethod("calc_metabolite_conservation", "metab_analyser", function(object, which_data, verbose=F) {
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
    
    id_split <- data.frame(
      met=rep(data_var,each=length(unique(data_merged$timepoint))),
      timepoint=rep(unique(data_merged$timepoint),length(data_var)),
      stringsAsFactors = F
    ) %>% 
      dplyr::mutate(id=paste0(met, "_",timepoint)) %>% 
      tidyr::spread(timepoint, id)
    
    
    if(verbose) cat("Processing dataframe", i,": ")
    out[[i]] = lapply(1:nrow(tp_split), function(x){
      
      if(verbose) cat(tp_split[x,] %>% paste0(collapse="vs"), "; ")
      
      
      v1_data <- data_merged %>% dplyr::filter(timepoint == tp_split[x,"V1"]) %>% 
        dplyr::select(-id, -timepoint,-rid)
      colnames(v1_data)  = paste0(names(v1_data), "_",tp_split$V1[x])
      
      v2_data <- data_merged %>% dplyr::filter(timepoint == tp_split[x,"V2"]) %>% 
        dplyr::select(-id, -timepoint,-rid) 
      colnames(v2_data)  = paste0(names(v2_data), "_",tp_split$V2[x])
      
      cor_mat = cbind(v1_data,v2_data) %>% 
        cor(use="pairwise.complete.obs")
      
      ci_out <- lapply(1:nrow(id_split), function(y){
        v2_id = id_split[,paste0(tp_split$V2[[x]])][y]
        v1_id = id_split[,paste0(tp_split$V1[[x]])][y]
        if(all(!is.na(c(v2_id, v1_id)))){
          ci_info<-cor_mat[intersect(rownames(cor_mat),id_split[,paste0(tp_split$V2[[x]])]),v1_id] %>% sort(decreasing = T)
          rank<-which(names(ci_info)==v2_id)
          
          data.frame(
            id=id_split$met[y],
            from_tp=tp_split[x,"V1"],
            to_tp=tp_split[x,"V2"],
            nsubject=length(ci_info),
            rank=ifelse(length(rank)>0, rank, NA),
            stringsAsFactors = F
          )
        }else{
          
        }
        
        
      }) %>% 
        do.call(what=rbind.data.frame)
      
      ci_out = ci_out %>%
        dplyr::mutate(ci=(1 - ((rank -1)/(nsubject-1)))) %>% 
        dplyr::arrange(ci) %>% 
        dplyr::mutate(y=ci, 
                      x=1:nrow(.[])) %>% 
        dplyr::select(x,y,ci,id, from_tp, to_tp, nsubject,rank)
      
      ci_out
    })
    
    names(out[[i]]) = paste0(tp_split$V1,"vs", tp_split$V2)
  }
  
  return(out)
})


#' Function to calculate dimensionality reduction methods such as tsne, umap and pca.
#' @description A method to apply on s4 object of class metab_analyse in order to obtain information after dimensionality reduction on a dataset/s
#' @param object An object of class metab_analyse
#' @param which_data a character vector - Names of the dataset from which the samples will be extracted
#' @param type type of the dimensionality reduction method to be applied. Accepted inputs are "UMAP", "tSNE", "PCA"
#' 
#' @return  a list with two dataframes containing the dimensionality reduction information 1) samples - data of the individuals(".$samples")
#'                     2) metabs - data of the metabolites(".$metabs")
setGeneric("calc_dimensionality_reduction", function(object, which_data, type) standardGeneric("calc_dimensionality_reduction"))

setMethod("calc_dimensionality_reduction", "metab_analyser", function(object, which_data, type) {
      if(length(which_data) > 1) {
        object@list_of_data <- mod_common_sample_extractor(object@list_of_data)
        data <- object@list_of_data[names(object@list_of_data) %in% which_data]
        data <- lapply(data, function(x) {
            return(x[sort(rownames(x)),])
          })
        data <- as.data.frame(do.call(cbind, data))
        colnames(data) <- unlist(lapply(strsplit(colnames(data), split="_data.", fixed=TRUE), function(x) return(x[2])))
      } else {
        data <- object@list_of_data[names(object@list_of_data) %in% which_data][[1]]
      }
      if(type %in% "PCA") {
        pca_metabs <- prcomp(t(data), scale.=T, center=T)
        pca_individuals <- prcomp(data, scale.=T, center=T)
        dr_data_metabs <- as.data.frame(pca_metabs$x[,1:2])
        dr_data_samples <- as.data.frame(pca_individuals$x[,1:2])
      } else if(type %in% "UMAP") {
        umap_individuals <- umap::umap(data)
        dr_data_samples <- as.data.frame(umap_individuals$layout)
        umap_metabs <- umap::umap(t(data))
        dr_data_metabs <- as.data.frame(umap_metabs$layout) 
      } else if(type %in% "tSNE") {
        tsne_samples <- tsne(t(data))
        tsne_metabs <- tsne(data)
        dr_data_metabs <- as.data.frame(tsne_metabs$data)
        dr_data_samples <- as.data.frame(tsne_samples$data)
      }
      return(list(metabs=dr_data_metabs, samples=dr_data_samples))
  })