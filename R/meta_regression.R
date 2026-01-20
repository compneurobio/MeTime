#' Meta comparison for regression outputs
#' @description Compare regression outputs within or across results.
#' @param object a S4 object of class metime_analyser or a list of two metime_analyser objects
#' @param method a character vector of methods 'sign', 'cor', 'het'
#' @param result_index character/numeric input for results. If NULL, all matching results are used.
#' @param name a character input to set the name of the results
#' @return An S4 object of class meta_results with the compared results and meta results
#' @export
setGeneric("meta_regression", function(object, method=c("sign", "cor", "het"), result_index=NULL, name="meta_regression_1") standardGeneric("meta_regression"))
meta_regression_impl <- function(object, method=c("sign", "cor", "het"), result_index=NULL, name="meta_regression_1") {
  method <- unique(method)
  if (!all(method %in% c("sign", "cor", "het"))) {
    stop('meta_regression(): method has to be one of "sign", "cor", or "het".')
  }
  analyzers <- meta_unpack_analyzers(object, function_name="meta_regression")
  results <- meta_collect_results(analyzers, result_index, allowed_calc_types="regression", function_name="meta_regression")
  comparisons <- meta_build_regression_comparisons(results)
  out <- list()
  for (i in seq_along(comparisons)) {
    comp_out <- meta_compare_regression(comparisons[[i]], method)
    if (is.null(comp_out)) {
      next
    }
    comp_names <- paste(names(comparisons)[i], names(comp_out), sep="__")
    out[comp_names] <- comp_out
  }
  return(meta_make_analyser(analyzers, results, out, calc_type="meta_regression",
                            calc_info=names(out), name=name, function_name="meta_regression",
                            params=list(result_index=result_index, method=method)))
}

setMethod("meta_regression", "metime_analyser", function(object, method=c("sign", "cor", "het"), result_index=NULL, name="meta_regression_1") {
  meta_regression_impl(object, method, result_index, name)
})

setMethod("meta_regression", "list", function(object, method=c("sign", "cor", "het"), result_index=NULL, name="meta_regression_1") {
  meta_regression_impl(object, method, result_index, name)
})