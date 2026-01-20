#' Meta comparison for matrix similarity
#' @description Compare pairwise distance or correlation results within or across results.
#' @param object a S4 object of class metime_analyser or a list of two metime_analyser objects
#' @param result_index character/numeric input for results. If NULL, all matching results are used.
#' @param name a character input to set the name of the results
#' @return An S4 object of class meta_results with the compared results and meta results
#' @export
setGeneric("meta_matrix_similarity", function(object, result_index=NULL, name="meta_matrix_similarity_1") standardGeneric("meta_matrix_similarity"))
meta_matrix_similarity_impl <- function(object, result_index=NULL, name="meta_matrix_similarity_1") {
  analyzers <- meta_unpack_analyzers(object, function_name="meta_matrix_similarity")
  results <- meta_collect_results(analyzers, result_index, allowed_calc_types=c("pairwise_distance", "pairwise_correlation"),
                                  function_name="meta_matrix_similarity")
  comparisons <- meta_get_comparison_builder()(results, compare_label="matrix_similarity")
  out <- list()
  for (i in seq_along(comparisons)) {
    comp <- comparisons[[i]]
    comp_out <- meta_compare_matrix_similarity(comp)
    if (!is.null(comp_out)) {
      out[[names(comparisons)[i]]] <- comp_out
    }
  }
  return(meta_make_analyser(analyzers, results, out, calc_type="meta_matrix_similarity",
                            calc_info=names(out), name=name, function_name="meta_matrix_similarity",
                            params=list(result_index=result_index)))
}

setMethod("meta_matrix_similarity", "metime_analyser", function(object, result_index=NULL, name="meta_matrix_similarity_1") {
  meta_matrix_similarity_impl(object, result_index, name)
})

setMethod("meta_matrix_similarity", "list", function(object, result_index=NULL, name="meta_matrix_similarity_1") {
  meta_matrix_similarity_impl(object, result_index, name)
})
