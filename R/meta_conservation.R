#' Meta comparison for conservation index results
#' @description Compare conservation index results within a result or across results.
#' @param object a S4 object of class metime_analyser or a list of two metime_analyser objects
#' @param result_index character/numeric input for results. If NULL, all matching results are used.
#' @param top_k numeric indicating the top-K features used for overlap calculations
#' @param name a character input to set the name of the results
#' @return An S4 object of class meta_results with the compared results and meta results
#' @export
setGeneric("meta_conservation", function(object, result_index=NULL, top_k=50, name="meta_conservation_1") standardGeneric("meta_conservation"))
meta_conservation_impl <- function(object, result_index=NULL, top_k=50, name="meta_conservation_1") {
  analyzers <- meta_unpack_analyzers(object, function_name="meta_conservation")
  results <- meta_collect_results(analyzers, result_index, allowed_calc_types=c("CI_metabolite", "CI_metabotype"),
                                  function_name="meta_conservation")
  comparisons <- meta_build_conservation_comparisons(results)
  out <- list(conservation_metabolite=list(), conservation_metabotype=list())
  plots <- list(conservation_metabolite=list(), conservation_metabotype=list())
  for (i in seq_along(comparisons)) {
    comp <- comparisons[[i]]
    comp_out <- meta_compare_conservation(comp, top_k=top_k)
    if (!is.null(comp_out)) {
      comp_name <- names(comparisons)[i]
      group_name <- meta_conservation_group_name(comp)
      for (metric in names(comp_out)) {
        if (is.null(out[[group_name]][[metric]])) {
          out[[group_name]][[metric]] <- list()
        }
        out[[group_name]][[metric]][[comp_name]] <- comp_out[[metric]]
      }
      plot_out <- meta_plot_conservation(comp, comp_out)
      if (length(plot_out) > 0) {
        for (metric in names(plot_out)) {
          if (is.null(plots[[group_name]][[metric]])) {
            plots[[group_name]][[metric]] <- list()
          }
          plots[[group_name]][[metric]][[comp_name]] <- plot_out[[metric]]
        }
      }
    }
  }
  out <- out[vapply(out, function(x) length(x) > 0, logical(1))]
  plots <- plots[names(out)]
  return(meta_make_analyser(analyzers, results, out, plots=plots, calc_type="meta_conservation",
                            calc_info=names(out), name=name, function_name="meta_conservation",
                            params=list(result_index=result_index, top_k=top_k)))
}

setMethod("meta_conservation", "metime_analyser", function(object, result_index=NULL, top_k=50, name="meta_conservation_1") {
  meta_conservation_impl(object, result_index, top_k, name)
})

setMethod("meta_conservation", "list", function(object, result_index=NULL, top_k=50, name="meta_conservation_1") {
  meta_conservation_impl(object, result_index, top_k, name)
})