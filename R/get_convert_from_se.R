#' Convert SummarizedExperiment objects back to metime_analyser
#'
#' This helper reverses the behaviour of `get_convert_to_se` by creating a new
#' `metime_analyser` from one or more `SummarizedExperiment` objects. It
#' reconstructs the original data matrix orientation (samples as rows and features as
#' columns) and retrieves row- and column-level metadata stored in the
#' `colData` and `rowData` of the SummarizedExperiment, respectively.
#' Datasets that lack the mandatory `id`, `subject` and `time` columns in
#' their row-level metadata will trigger an error, as these fields are required
#' for a valid `metime_analyser`.
#'
#' @param se A single `SummarizedExperiment` or a list of such objects. When a list
#'   is provided, each entry will be added as a separate dataset to the returned
#'   `metime_analyser` object.
#' @param name Optional character vector of dataset names. If `se` is a list and
#'   `name` is `NULL`, the names of the list elements will be used. When
#'   `se` is a single object and `name` is `NULL`, the default
#'   dataset name "set_1" will be applied.
#' @param annotations_index A named list specifying phenotype and medication datasets.
#'   This mirrors the argument in `get_make_analyser_object`. When converting
#'   multiple datasets, the same annotations list will be applied to all.
#' @return An object of class `metime_analyser` containing one or more
#'   reconstructed datasets.
#' @examples
#' \dontrun{
#'   # Assume 'se_obj' is a SummarizedExperiment produced by get_convert_to_se()
#'   analyser <- get_convert_from_se(se_obj, name = "my_dataset")
#' }
#' @export
setGeneric("get_convert_from_se", function(se, name=NULL, annotations_index=list()) {
  standardGeneric("get_convert_from_se")
})

## Method for a single SummarizedExperiment object
setMethod("get_convert_from_se", "SummarizedExperiment",
          function(se, name=NULL, annotations_index=list()) {
            if(is.null(name)) {
              name <- "set_1"
            }
            # extract assay, rowData (feature metadata) and colData (sample metadata)
            assay_data <- SummarizedExperiment::assay(se)
            row_meta <- as.data.frame(SummarizedExperiment::rowData(se))
            col_meta <- as.data.frame(SummarizedExperiment::colData(se))
            # Transpose assay back to samples x features
            data_mat <- t(as.matrix(assay_data))
            # ensure that metadata contains id column matching row/col names
            if(!("id" %in% names(col_meta))) {
              stop("Column-level metadata must contain an 'id' column matching column names of the SummarizedExperiment.")
            }
            if(!("id" %in% names(row_meta))) {
              stop("Row-level metadata must contain an 'id' column matching row names of the SummarizedExperiment.")
            }
            # reorder metadata to match data matrix
            col_meta <- col_meta[match(colnames(assay_data), col_meta$id), , drop=FALSE]
            row_meta <- row_meta[match(rownames(assay_data), row_meta$id), , drop=FALSE]
            # create metime analyser; validity requires subject and time fields
            required_fields <- c("id","subject","time")
            if(!all(required_fields %in% names(col_meta))) {
              stop("The column-level metadata must contain 'id', 'subject' and 'time' columns for a valid metime_analyser.")
            }
            # Build the analyser object
            out <- get_make_analyser_object(data=data_mat,
                                            col_data=row_meta,
                                            row_data=col_meta,
                                            annotations_index=annotations_index,
                                            name=name,
                                            results=list())
            return(out)
          })

## Method for list of SummarizedExperiment objects
#' @describeIn get_convert_from_se Convert a list of SummarizedExperiment objects.
#' When a single dataset has been split by timepoint (i.e. the list contains
#' multiple SummarizedExperiments whose names share a common prefix before the
#' final underscore), this method automatically combines them back into a
#' single dataset by row-binding their sample-level data and taking the
#' feature-level metadata from the first element.
setMethod("get_convert_from_se", "list",
          function(se, name=NULL, annotations_index=list()) {
            if(length(se) == 0) stop("The list of SummarizedExperiment objects is empty.")
            # ensure names are present
            if(is.null(names(se))) {
              stop("The list of SummarizedExperiment objects must have names to infer dataset grouping.")
            }
            dataset_names <- names(se)
            # group datasets by base name (prefix before the last underscore)
            base_names <- sub("_([^_]+)$", "", dataset_names)
            unique_bases <- unique(base_names)
            # prepare list to hold combined datasets
            combined <- list()
            for(bn in unique_bases) {
              indices <- which(base_names == bn)
              if(length(indices) == 1) {
                # single SE; convert directly
                se_obj <- se[[indices]]
                combined[[bn]] <- get_convert_from_se(se_obj, name=bn, annotations_index=annotations_index)
              } else {
                # multiple SEs for the same base; combine
                data_mats <- list()
                row_meta_list <- list()
                col_meta <- NULL
                for(idx in indices) {
                  se_i <- se[[idx]]
                  assay_data <- SummarizedExperiment::assay(se_i)
                  row_meta <- as.data.frame(SummarizedExperiment::colData(se_i))
                  col_meta_i <- as.data.frame(SummarizedExperiment::rowData(se_i))
                  # transpose assay to sample x feature
                  data_mats[[length(data_mats)+1]] <- t(as.matrix(assay_data))
                  row_meta_list[[length(row_meta_list)+1]] <- row_meta
                  if(is.null(col_meta)) col_meta <- col_meta_i
                }
                # combine sample-level data and metadata
                combined_data <- do.call(rbind, data_mats)
                combined_row_meta <- dplyr::bind_rows(row_meta_list)
                # re-order feature metadata to match combined_data columns
                col_meta <- col_meta[match(colnames(combined_data), col_meta$id), , drop=FALSE]
                combined_row_meta <- combined_row_meta[match(rownames(combined_data), combined_row_meta$id), , drop=FALSE]
                # Build analyser for this combined dataset
                analyser_tmp <- get_make_analyser_object(data=combined_data,
                                                         col_data=col_meta,
                                                         row_data=combined_row_meta,
                                                         annotations_index=annotations_index,
                                                         name=bn,
                                                         results=list())
                combined[[bn]] <- analyser_tmp
              }
            }
            # merge combined analysers into a single analyser if necessary
            # start with first element
            analyser <- combined[[1]]
            if(length(combined) > 1) {
              cnames <- names(combined)
              # skip first
              for(j in 2:length(combined)) {
                an_tmp <- combined[[j]]
                ds_name <- names(an_tmp@list_of_data)[1]
                analyser <- add_dataset(object=analyser,
                                        data=an_tmp@list_of_data[[ds_name]],
                                        col_data=an_tmp@list_of_col_data[[ds_name]],
                                        row_data=an_tmp@list_of_row_data[[ds_name]],
                                        name=ds_name)
              }
            }
            return(analyser)
          })

## Method for MultiAssayExperiment objects
#' @rdname get_convert_from_se
#' @export
setMethod("get_convert_from_se", "MultiAssayExperiment",
          function(se, name=NULL, annotations_index=list()) {
            if(!requireNamespace("MultiAssayExperiment", quietly=TRUE)) {
              stop("MultiAssayExperiment package is required to convert from MultiAssayExperiment objects.")
            }
            # extract experiments
            exps <- MultiAssayExperiment::experiments(se)
            if(length(exps) == 0) stop("The MultiAssayExperiment contains no experiments to convert.")
            # determine dataset names
            dataset_names <- names(exps)
            if(is.null(name)) {
              name <- dataset_names
            }
            if(length(name) != length(exps)) {
              stop("Please provide names for each experiment in the MultiAssayExperiment.")
            }
            # convert experiments individually using list method
            se_list <- as.list(exps)
            names(se_list) <- name
            return(get_convert_from_se(se_list, name=name, annotations_index=annotations_index))
          })