#' Convert analyser datasets to SummarizedExperiment objects
#' @description Converts each dataset in a metime_analyser into a SummarizedExperiment object.
#' Phenotype and medication datasets are skipped. Longitudinal datasets are split into
#' separate SummarizedExperiment objects by timepoint with a timepoint suffix.
#' @param object a S4 object of class "metime_analyser".
#' @param which_data character vector of dataset names to convert. If NULL, converts all datasets.
#' @param exclude character vector of dataset names to skip. Defaults to phenotype/medication data.
#' @return A named list of SummarizedExperiment objects.
#' @export
setGeneric("get_convert_to_se", function(object, which_data=NULL, exclude=c("phenotype_data", "medication_data")) {
  standardGeneric("get_convert_to_se")
})

setMethod("get_convert_to_se", "metime_analyser", function(object, which_data=NULL, exclude=c("phenotype_data", "medication_data")) {
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

    timepoints <- if("time" %in% names(row_data)) unique(row_data$time) else NA
    timepoints <- timepoints[!is.na(timepoints)]
    if(length(timepoints) == 0) timepoints <- NA

    if(length(timepoints) > 1) {
      for(timepoint in timepoints) {
        keep <- row_data$time %in% timepoint
        subset_data <- data[keep, , drop=FALSE]
        subset_row_data <- row_data[keep, , drop=FALSE]
        se_name <- paste0(dataset, "_", as.character(timepoint))
        out[[se_name]] <- build_se_object(subset_data, subset_row_data, col_data)
      }
    } else {
      out[[dataset]] <- build_se_object(data, row_data, col_data)
    }
  }

  out
})

build_se_object <- function(data, row_data, col_data) {
  assay <- t(as.matrix(data))
  if("id" %in% names(row_data)) {
    row_data <- row_data[match(rownames(data), row_data$id), , drop=FALSE]
  } else {
    row_data <- row_data[rownames(data), , drop=FALSE]
  }
  if("id" %in% names(col_data)) {
    col_data <- col_data[match(colnames(data), col_data$id), , drop=FALSE]
  } else {
    col_data <- col_data[colnames(data), , drop=FALSE]
  }
  rownames(row_data) <- rownames(data)
  rownames(col_data) <- colnames(data)

  SummarizedExperiment::SummarizedExperiment(
    assays = list(data = assay),
    rowData = S4Vectors::DataFrame(col_data),
    colData = S4Vectors::DataFrame(row_data)
  )
}
