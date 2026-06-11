#' Convert a metime_analyser into a MultiAssayExperiment
#'
#' This function constructs a `MultiAssayExperiment` from the datasets stored in a
#' `metime_analyser`. Each dataset is first converted into one or more
#' `SummarizedExperiment` objects via the `get_convert_to_se` helper. When
#' `split_by_time` is `TRUE`, longitudinal datasets are split by their
#' timepoints into separate experiments; otherwise, the full dataset is
#' preserved in a single `SummarizedExperiment` with a time column in the
#' `colData`. Phenotype and medication datasets are skipped by default.
#'
#' @param object A `metime_analyser` object containing one or more datasets.
#' @param which_data Character vector of dataset names to convert. If `NULL`,
#'   all datasets are converted. Data sets listed in `exclude` will be
#'   ignored.
#' @param exclude Character vector of dataset names to skip. Defaults to
#'   `c("phenotype_data", "medication_data")`.
#' @param split_by_time Logical indicating whether longitudinal datasets
#'   should be split by timepoint. Passed through to `get_convert_to_se`.
#'   Default is `TRUE`.
#' @return A `MultiAssayExperiment` containing one `SummarizedExperiment` per
#'   selected dataset (or timepoint-specific split). Sample-level metadata is
#'   not merged across datasets, mirroring the structure of the original
#'   analyser object.
#' @examples
#' \dontrun{
#'   # convert to MAE with splitting by time (default)
#'   mae <- get_convert_to_mae(humet_object)
#'   mae
#'
#'   # convert to MAE without splitting by time
#'   mae_full <- get_convert_to_mae(humet_object, split_by_time = FALSE)
#'   mae_full
#' }
#' @export
setGeneric("get_convert_to_mae", function(object, which_data=NULL,
                                          exclude=c("phenotype_data", "medication_data"),
                                          split_by_time=TRUE) {
  standardGeneric("get_convert_to_mae")
})

setMethod("get_convert_to_mae", "metime_analyser",
          function(object, which_data=NULL,
                   exclude=c("phenotype_data", "medication_data"),
                   split_by_time=TRUE) {
            if(!requireNamespace("MultiAssayExperiment", quietly=TRUE)) {
              stop("MultiAssayExperiment package is required to convert to a MultiAssayExperiment.")
            }
            # convert datasets to SummarizedExperiment objects using helper; honour split_by_time
            se_list <- get_convert_to_se(object,
                                         which_data=which_data,
                                         exclude=exclude,
                                         split_by_time=split_by_time)
            # wrap into an ExperimentList; SummarizedExperiment list may already carry
            # shared rowData objects across timepoints to avoid duplication
            exp_list <- MultiAssayExperiment::ExperimentList(se_list)
            # create MAE; sampleMap and colData will be inferred from experiment colnames
            mae <- MultiAssayExperiment::MultiAssayExperiment(experiments=exp_list)
            return(mae)
          })