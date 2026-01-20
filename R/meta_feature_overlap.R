#' Meta comparison for feature overlap
#' @description Compare feature selection outputs within or across results.
#' @param object a S4 object of class metime_analyser or a list of two metime_analyser objects
#' @param result_index character/numeric input for results. If NULL, all matching results are used.
#' @param name a character input to set the name of the results
#' @return An S4 object of class meta_results with the compared results and meta results
#' @export
setGeneric("meta_feature_overlap", function(object, result_index=NULL, name="meta_feature_overlap_1") standardGeneric("meta_feature_overlap"))
meta_feature_overlap_impl <- function(object, result_index=NULL, name="meta_feature_overlap_1") {
  analyzers <- meta_unpack_analyzers(object, function_name="meta_feature_overlap")
  results <- meta_collect_results(analyzers, result_index, allowed_calc_types="feature_selection", function_name="meta_feature_overlap")
  comparisons <- meta_get_comparison_builder()(results, compare_label="feature_overlap")
  out <- lapply(seq_along(comparisons), function(i) meta_compare_feature_overlap(comparisons[[i]]))
  names(out) <- names(comparisons)
  out <- out[!vapply(out, is.null, logical(1))]
  return(meta_make_analyser(analyzers, results, out, calc_type="meta_feature_overlap",
                            calc_info=names(out), name=name, function_name="meta_feature_overlap",
                            params=list(result_index=result_index)))
}

setMethod("meta_feature_overlap", "metime_analyser", function(object, result_index=NULL, name="meta_feature_overlap_1") {
  meta_feature_overlap_impl(object, result_index, name)
})

setMethod("meta_feature_overlap", "list", function(object, result_index=NULL, name="meta_feature_overlap_1") {
  meta_feature_overlap_impl(object, result_index, name)
})