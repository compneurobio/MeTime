#' Convert a MultiAssayExperiment into a metime_analyser
#'
#' This is a convenience wrapper around `get_convert_from_se` that converts a
#' `MultiAssayExperiment` (MAE) into a `metime_analyser`. Each experiment in the
#' MAE is first converted into a `SummarizedExperiment` and then combined into
#' a single `metime_analyser` with one dataset per experiment. If dataset
#' names are not provided, the names of the experiments are used.
#'
#' @param mae A `MultiAssayExperiment` object.
#' @param name Optional character vector of dataset names corresponding to the
#'   experiments within the MAE. If `NULL`, the experiment names will be used.
#' @param annotations_index A named list specifying phenotype and medication
#'   datasets, passed through to the underlying conversion routines.
#' @return A `metime_analyser` object containing one dataset per experiment in
#'   the input MAE.
#' @examples
#' \dontrun{
#'   mae <- get_convert_to_mae(humet_object)
#'   analyser <- get_convert_from_mae(mae)
#' }
#' @export
setGeneric("get_convert_from_mae", function(mae, name=NULL, annotations_index=list()) {
  standardGeneric("get_convert_from_mae")
})

setMethod("get_convert_from_mae", "MultiAssayExperiment",
          function(mae, name=NULL, annotations_index=list()) {
            # delegate to the existing get_convert_from_se method for MAE
            get_convert_from_se(mae, name=name, annotations_index=annotations_index)
          })