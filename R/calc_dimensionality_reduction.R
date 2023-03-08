
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
#' @param cols_for_meta A Character vector to define columns names that are to be used for plotting purposes
#' @param stratifications List to stratify data into a subset. Usage list(name=value)
#' @param ... additional arguments that can be passed on to prcomp(), M3C::tsne() and umap::umap()
#' @importClassesFrom metime_analyser
#' @return  a list with two plotter objects containing the dimensionality reduction information that can be parsed into plotting function 
#'                   
#' @export
setGeneric("calc_dimensionality_reduction_samples", function(object, which_data, type, cols_for_meta, stratifications=NULL, ...) standardGeneric("calc_dimensionality_reduction_samples"))

setMethod("calc_dimensionality_reduction_samples", "metime_analyser", function(object, which_data, type, cols_for_meta, stratifications=NULL, ...) {
      
      data_list <- get_stratified_data(object=object, which_data=which_data, stratifications=stratifications)
      data <- data_list[["data"]]
      row_data <- data_list[["row_data"]]
      metadata_samples <- get_metadata_for_rows(object=object, which_data=which_data, columns=cols_for_samples)
      if(type %in% "PCA") {
        pca_individuals <- prcomp(data, retx=TRUE, scale.=TRUE, center=TRUE, tol = NULL, ...)
        dr_data_samples <- as.data.frame(pca_individuals$x[,1:2])
        out <- get_make_results(object=object, data=list(dr_data_samples), 
                metadata=list(metabs=metadata_metabs, samples=metadata_samples),
                calc_type="PCA_samples", 
                calc_info=paste("dimensionality reduction method:", 
                    type, "for samples data of", 
                    paste(which_data, collapse=" & "), sep=" "),
                name=paste("PCA_samples", paste(which_data, collapse="_&_"), sep=""))
      } else if(type %in% "UMAP") {
        umap_individuals <- umap::umap(data, ...)
        dr_data_samples <- as.data.frame(umap_individuals$layout)
        colnames(dr_data_samples) <- c("UMAP1", "UMAP2")
        out <- get_make_results(object=object, data=list(dr_data_samples), 
                metadata=list(metabs=metadata_metabs, samples=metadata_samples),
                calc_type="UMAP_samples", 
                calc_info=paste("dimensionality reduction method:", 
                    type, "for samples data of", 
                    paste(which_data, collapse=" & "), sep=" "),
                name=paste("UMAP_samples", paste(which_data, collapse="_&_"), sep=""))
      } else if(type %in% "tSNE") {
        tsne_samples <- M3C::tsne(t(data), ...)
        dr_data_samples <- as.data.frame(tsne_samples$data)
        rownames(dr_data_samples) <- rownames(data)
        out <- get_make_results(object=object, data=list(dr_data_samples), 
                metadata=list(metabs=metadata_metabs, samples=metadata_samples),
                calc_type="tSNE_samples", 
                calc_info=paste("dimensionality reduction method:", 
                    type, "for samples data of", 
                    paste(which_data, collapse=" & "), sep=" "),
                name=paste("tSNE_samples", paste(which_data, collapse="_&_"), sep=""))
      }
      out <- add_function_info(object=out, function_name="calc_dimensionality_reduction_samples", 
        params=list(which_data=which_data, type=type, cols_for_meta=cols_for_meta, 
            ...))
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
#' @param cols_for_meta A list of Character vectors to define columns names that are to be used for plotting purposes
#' @param stratifications List to stratify data into a subset. Usage list(name=value)
#' @param ... additional arguments that can be passed on to prcomp(), M3C::tsne() and umap::umap()
#' @return  a list with two plotter objects containing the dimensionality reduction information that can be parsed into plotting function 
#'                   
#' @export
setGeneric("calc_dimensionality_reduction_metabs", function(object, which_data, type, cols_for_meta, stratifications=list(), ...) standardGeneric("calc_dimensionality_reduction_metabs"))

setMethod("calc_dimensionality_reduction_metabs", "metime_analyser", function(object, which_data, type, cols_for_meta, stratifications=list(), ...) {
      
      data_list <- get_stratified_data(object=object, which_data=which_data, stratifications=stratifications)
      data <- data_list[["data"]]
      row_data <- data_list[["row_data"]]
      metadata_metabs <- get_metadata_for_columns(object=object, 
        which_data=which_data, columns=cols_for_metabs)
      if(type %in% "PCA") {
        pca_metabs <- prcomp(t(data), retx=TRUE, scale.=TRUE, center=TRUE, tol = NULL, ...)
        dr_data_metabs <- as.data.frame(pca_metabs$x[,1:2])
        out <- get_make_results(object=object, data=list(dr_data_metabs), 
                metadata=metadata_metabs,
                calc_type="PCA_metabs", 
                calc_info=paste("dimensionality reduction method:", 
                    type, "for metabolites data of", 
                    paste(which_data, collapse=" & "), sep=" "),
                name=paste("PCA_metabs", paste(which_data, collapse="_&_"), sep=""))
      } else if(type %in% "UMAP") {
        umap_metabs <- umap::umap(t(data), ...)
        dr_data_metabs <- as.data.frame(umap_metabs$layout)
        colnames(dr_data_metabs) <- c("UMAP1", "UMAP2")
        out <- get_make_results(object=object, data=list(dr_data_metabs), 
                metadata=metadata_metabs,
                calc_type="UMAP_metabs", 
                calc_info=paste("dimensionality reduction method:", 
                    type, "for metabolites data of", 
                    paste(which_data, collapse=" & "), sep=" "),
                name=paste("UMAP_metabs", paste(which_data, collapse="_&_"), sep=""))
      } else if(type %in% "tSNE") {
        tsne_metabs <- M3C::tsne(data, ...)
        dr_data_metabs <- as.data.frame(tsne_metabs$data)
        rownames(dr_data_metabs) <- colnames(data)
        out <- get_make_results(object=object, data=list(dr_data_metabs), 
                metadata=metadata_metabs,
                calc_type="tSNE", 
                calc_info=paste("dimensionality reduction method:", 
                    type, "for metabolites data of", 
                    paste(which_data, collapse=" & "), sep=" "),
                name=paste("tSNE_", paste(which_data, collapse="_&_"), sep=""))
      }
      out <- add_function_info(object=out, function_name="calc_dimensionality_reduction_metabs", 
        params=list(which_data=which_data, type=type, cols_for_meta=cols_for_meta, 
            ...))
      return(out)
  })
