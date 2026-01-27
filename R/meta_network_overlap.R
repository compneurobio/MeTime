

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
  out <- list(network_overlap=list(jaccard=list()))
  plots <- list(network_overlap=list(jaccard=list()))
  for (i in seq_along(comparisons)) {
    comp_out <- meta_compare_network(comparisons[[i]])
    if (!is.null(comp_out)) {
      comp_name <- names(comparisons)[i]
      out$network_overlap$jaccard[[comp_name]] <- comp_out
      plots$network_overlap$jaccard[[comp_name]] <- meta_plot_network_overlap(comparisons[[i]])
    }
  }
  out <- out[vapply(out, function(x) length(x) > 0, logical(1))]
  plots <- plots[names(out)]
  return(meta_make_analyser(analyzers, results, out, plots=plots, calc_type="meta_network_overlap",
                            calc_info=names(out), name=name, function_name="meta_network_overlap",
                            params=list(result_index=result_index)))
}

setMethod("meta_network_overlap", "metime_analyser", function(object, result_index=NULL, name="meta_network_overlap_1") {
  meta_network_overlap_impl(object, result_index, name)
})

setMethod("meta_network_overlap", "list", function(object, result_index=NULL, name="meta_network_overlap_1") {
  meta_network_overlap_impl(object, result_index, name)
})