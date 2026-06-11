#' Convert analyser datasets to SummarizedExperiment objects with optional time splitting
#'
#' @description
#' Converts each dataset in a `metime_analyser` into one or more
#' `SummarizedExperiment` objects. When `split_by_time` is `TRUE`, longitudinal
#' datasets (those whose row-level metadata include a `time` column with more
#' than one unique value) are split into separate `SummarizedExperiment`s for
#' each timepoint. Otherwise, the entire dataset is stored in a single
#' `SummarizedExperiment` with the `time` variable preserved in the
#' `colData` slot. Feature metadata (`rowData`) is not duplicated across
#' timepoints; the same `DataFrame` is reused for each split.
#'
#' Phenotype and medication datasets are skipped by default.
#'
#' @param object A `metime_analyser` object.
#' @param which_data Character vector of dataset names to convert. If `NULL`,
#'   all datasets are converted.
#' @param exclude Character vector of dataset names to skip. Defaults to
#'   `c("phenotype_data", "medication_data")`.
#' @param split_by_time Logical indicating whether longitudinal datasets should
#'   be split by timepoint. If `FALSE`, a single `SummarizedExperiment` is
#'   returned per dataset with a `time` column in the sample-level metadata.
#' @return A named list of `SummarizedExperiment` objects.
#' @export
setGeneric("get_convert_to_se", function(object, which_data=NULL,
                                          exclude=c("phenotype_data", "medication_data"),
                                          split_by_time=TRUE) {
  standardGeneric("get_convert_to_se")
})

setMethod("get_convert_to_se", "metime_analyser",
          function(object, which_data=NULL, exclude=c("phenotype_data", "medication_data"),
                   split_by_time=TRUE) {
            if(!requireNamespace("SummarizedExperiment", quietly=TRUE)) {
              stop("SummarizedExperiment package is required to convert data to SE objects.")
            }
            if(!requireNamespace("S4Vectors", quietly=TRUE)) {
              stop("S4Vectors package is required to convert data to SE objects.")
            }
            data_names <- names(object@list_of_data)
            if(is.null(which_data)) which_data <- data_names
            missing <- setdiff(which_data, data_names)
            if(length(missing) > 0) warning("Datasets not found: ", paste(missing, collapse=", "))
            which_data <- intersect(which_data, data_names)
            out <- list()
            for(dataset in which_data) {
              if(dataset %in% exclude) next
              data <- object@list_of_data[[dataset]]
              row_data <- object@list_of_row_data[[dataset]]
              col_data <- object@list_of_col_data[[dataset]]
              if(is.null(data) || is.null(row_data) || is.null(col_data)) {
                warning("Missing data, row_data, or col_data for dataset: ", dataset)
                next
              }
              # Determine timepoints from sample-level metadata (row_data)
              timepoints <- if("time" %in% names(row_data)) unique(row_data$time) else NA
              timepoints <- timepoints[!is.na(timepoints)]
              if(length(timepoints) == 0) timepoints <- NA
              # feature names and sample names for ordering
              feature_names <- colnames(data)
              sample_names <- rownames(data)
              # order and label feature-level metadata to match assay rows (features)
              # reorder by 'id' if present
              col_data_ord <- col_data
              if("id" %in% names(col_data)) {
                col_data_ord <- col_data[match(feature_names, col_data$id), , drop=FALSE]
              } else if(!is.null(rownames(col_data))) {
                col_data_ord <- col_data[match(feature_names, rownames(col_data)), , drop=FALSE]
              }
              rownames(col_data_ord) <- feature_names
              feature_df <- S4Vectors::DataFrame(col_data_ord)
              # If split_by_time and multiple timepoints, create separate SE objects per timepoint
              if(split_by_time && length(timepoints) > 1) {
                for(tp in timepoints) {
                  keep <- row_data$time %in% tp
                  subset_data <- data[keep, , drop=FALSE]
                  subset_row_data <- row_data[keep, , drop=FALSE]
                  # sample names for this timepoint
                  sample_tp <- rownames(subset_data)
                  # reorder sample metadata to match assay columns
                  row_data_ord <- subset_row_data
                  if("id" %in% names(subset_row_data)) {
                    row_data_ord <- subset_row_data[match(sample_tp, subset_row_data$id), , drop=FALSE]
                  } else if(!is.null(rownames(subset_row_data))) {
                    row_data_ord <- subset_row_data[match(sample_tp, rownames(subset_row_data)), , drop=FALSE]
                  }
                  rownames(row_data_ord) <- sample_tp
                  sample_df <- S4Vectors::DataFrame(row_data_ord)
                  # build assay matrix (features x samples)
                  assay <- t(as.matrix(subset_data))
                  se_name <- paste0(dataset, "_", as.character(tp))
                  out[[se_name]] <- SummarizedExperiment::SummarizedExperiment(
                    assays = list(data = assay),
                    rowData = feature_df,
                    colData = sample_df
                  )
                }
              } else {
                # keep full dataset; include time variable in sample metadata
                # reorder row_data to match sample_names
                row_data_ord <- row_data
                if("id" %in% names(row_data)) {
                  row_data_ord <- row_data[match(sample_names, row_data$id), , drop=FALSE]
                } else if(!is.null(rownames(row_data))) {
                  row_data_ord <- row_data[match(sample_names, rownames(row_data)), , drop=FALSE]
                }
                rownames(row_data_ord) <- sample_names
                sample_df <- S4Vectors::DataFrame(row_data_ord)
                assay <- t(as.matrix(data))
                se_name <- dataset
                out[[se_name]] <- SummarizedExperiment::SummarizedExperiment(
                  assays = list(data = assay),
                  rowData = feature_df,
                  colData = sample_df
                )
              }
            }
            return(out)
          })