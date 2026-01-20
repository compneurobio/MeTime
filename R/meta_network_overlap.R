#' Meta comparison for network overlap
#' @description Compare network edge overlap within or across results.
#' @param object a S4 object of class metime_analyser or a list of two metime_analyser objects
#' @param result_index character/numeric input for results. If NULL, all matching results are used.
#' @param name a character input to set the name of the results
#' @return An S4 object of class meta_results with the compared results and meta results
#' @export
setGeneric("meta_network_overlap", function(object, result_index=NULL, name="meta_network_overlap_1") standardGeneric("meta_network_overlap"))
meta_network_overlap_impl <- function(object, result_index=NULL, name="meta_network_overlap_1") {
  analyzers <- meta_unpack_analyzers(object, function_name="meta_network_overlap")
  results <- meta_collect_results(analyzers, result_index, allowed_calc_types=c("genenet_ggm", "multibipartite_ggm", "temporal_network"),
                                  function_name="meta_network_overlap")
  comparisons <- meta_get_comparison_builder()(results, compare_label="network_overlap", allow_network_mismatch=TRUE)
  out <- lapply(seq_along(comparisons), function(i) meta_compare_network(comparisons[[i]]))
  names(out) <- names(comparisons)
  out <- out[!vapply(out, is.null, logical(1))]
  return(meta_make_analyser(analyzers, results, out, calc_type="meta_network_overlap",
                            calc_info=names(out), name=name, function_name="meta_network_overlap",
                            params=list(result_index=result_index)))
}

setMethod("meta_network_overlap", "metime_analyser", function(object, result_index=NULL, name="meta_network_overlap_1") {
  meta_network_overlap_impl(object, result_index, name)
})

setMethod("meta_network_overlap", "list", function(object, result_index=NULL, name="meta_network_overlap_1") {
  meta_network_overlap_impl(object, result_index, name)
})
