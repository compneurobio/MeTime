#' Meta comparison for conservation index results
#' @description Compare conservation index results within a result or across results.
#' @param object a S4 object of class metime_analyser or a list of two metime_analyser objects
#' @param result_index character/numeric input for results. If NULL, all matching results are used.
#' @param top_k numeric indicating the top-K features used for overlap calculations
#' @param name a character input to set the name of the results
#' @return An S4 object of class meta_analyser with the compared results and meta results
#' @export
setGeneric("meta_conservation", function(object, result_index=NULL, top_k=50, name="meta_conservation_1") standardGeneric("meta_conservation"))
setMethod("meta_conservation", "metime_analyser", function(object, result_index=NULL, top_k=50, name="meta_conservation_1") {
     analyzers <- meta_unpack_analyzers(object, function_name="meta_conservation")
  results <- meta_collect_results(analyzers, result_index, allowed_calc_types=c("CI_metabolite", "CI_metabotype"),
                                  function_name="meta_conservation")
  comparisons <- meta_build_conservation_comparisons(results)
  out <- list()
  for (i in seq_along(comparisons)) {
    comp_out <- meta_compare_conservation(comparisons[[i]], top_k=top_k)
    comp_names <- paste(names(comparisons)[i], names(comp_out), sep="__")
    out[comp_names] <- comp_out
  }
  names(out) <- make.unique(names(out))
  return(meta_make_analyser(analyzers, results, out, calc_type="meta_conservation",
                            calc_info=names(out), name=name, function_name="meta_conservation",
                            params=list(result_index=result_index, top_k=top_k)))
})

