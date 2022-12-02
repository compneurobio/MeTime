
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
#' @param cols_for_metabs a list of character vectors for getting metadata for columns for plotting purposes
#' @param cols_for_samples a character vector to define the columns to extract metadata for plotting purposes  
#' @param ... additional arguments that can be passed on to prcomp(), M3C::tsne() and umap::umap()
#' @return  a list with two plotter objects containing the dimensionality reduction information that can be parsed into plotting function 
#'                     1) samples - data of the individuals(".$samples")
#'                     2) metabs - data of the metabolites(".$metabs")
#' @export
setGeneric("calc_dimensionality_reduction", function(object, which_data, type, cols_for_metabs, cols_for_samples, ...) standardGeneric("calc_dimensionality_reduction"))

setMethod("calc_dimensionality_reduction", "metime_analyser", function(object, which_data, type, cols_for_metabs, cols_for_samples, ...) {
      if(length(which_data) > 1) {
        object <- mod_extract_common_samples(object)
        data <- object@list_of_data[names(object@list_of_data) %in% which_data]
        data <- lapply(data, function(x) {
            return(x[sort(rownames(x)),])
          })
        data <- as.data.frame(do.call(cbind, unname(data)))
      } else {
        data <- object@list_of_data[names(object@list_of_data) %in% which_data][[1]]
      }
      rm_col = intersect(names(data), c("time","subject"))
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
        tsne_samples <- M3C::tsne(t(data), ...)
        tsne_metabs <- M3C::tsne(data, ...)
        dr_data_metabs <- as.data.frame(tsne_metabs$data)
        rownames(dr_data_metabs) <- colnames(data)
        dr_data_samples <- as.data.frame(tsne_samples$data)
        rownames(dr_data_samples) <- rownames(data)
      }
      metadata_metabs <- get_metadata_for_columns(object=object, which_data=which_data, columns=cols_for_metabs, 
                 names=c("name", "pathway"), index_of_names=rep("id", each=length(which_data)))
      metadata_samples <- get_metadata_for_rows(object=object, which_data=which_data, columns=cols_for_samples)
      plotter_metabs <- get_make_plotter_object(data=dr_data_metabs, metadata=metadata_metabs,
                calc_type="Dimensionality Reduction", 
                calc_info=paste("dimensionality reduction method:", type, "for metabolites data", paste(which_data, collapse=" & "), sep=" "), 
                plot_type="dot", style="ggplot")
      plotter_samples <- get_make_plotter_object(data=dr_data_samples, metadata=metadata_samples,
                calc_type="Dimensionality Reduction", 
                calc_info=paste("dimensionality reduction method:", type, "for samples data", paste(which_data, collapse=" & "), sep=" "), 
                plot_type="dot", style="ggplot")
      out <- list(metabs=plotter_metabs, samples=plotter_samples)
      return(out)
  })


