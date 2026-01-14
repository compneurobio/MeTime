#' Meta comparison for conservation index results
#' @description Compare conservation index results within a result or across results.
#' @param object a S4 object of class metime_analyser or a list of two metime_analyser objects
#' @param result_index character/numeric input for results. If NULL, all matching results are used.
#' @param name a character input to set the name of the results
#' @return An S4 object of class meta_analyser with the compared results and meta results
#' @export
setGeneric("meta_conservation", function(object, result_index=NULL, name="meta_conservation_1") standardGeneric("meta_conservation"))
meta_conservation_impl <- function(object, result_index=NULL, name="meta_conservation_1") {
  analyzers <- meta_unpack_analyzers(object, function_name="meta_conservation")
  results <- meta_collect_results(analyzers, result_index, allowed_calc_types=c("CI_metabolite", "CI_metabotype"),
                                  function_name="meta_conservation")
  comparisons <- meta_build_conservation_comparisons(results)
  out <- lapply(seq_along(comparisons), function(i) meta_compare_conservation(comparisons[[i]])) %>%
    setNames(names(comparisons))
  return(meta_make_analyser(analyzers, results, out, calc_type="meta_conservation",
                            calc_info=names(out), name=name, function_name="meta_conservation",
                            params=list(result_index=result_index)))
}

setMethod("meta_conservation", "metime_analyser", function(object, result_index=NULL, name="meta_conservation_1") {
  meta_conservation_impl(object, result_index, name)
})

setMethod("meta_conservation", "list", function(object, result_index=NULL, name="meta_conservation_1") {
  meta_conservation_impl(object, result_index, name)
})

#' Meta comparison for matrix similarity
#' @description Compare pairwise distance or correlation results within or across results.
#' @param object a S4 object of class metime_analyser or a list of two metime_analyser objects
#' @param result_index character/numeric input for results. If NULL, all matching results are used.
#' @param name a character input to set the name of the results
#' @return An S4 object of class meta_analyser with the compared results and meta results
#' @export
setGeneric("meta_matrix_similarity", function(object, result_index=NULL, name="meta_matrix_similarity_1") standardGeneric("meta_matrix_similarity"))
meta_matrix_similarity_impl <- function(object, result_index=NULL, name="meta_matrix_similarity_1") {
  analyzers <- meta_unpack_analyzers(object, function_name="meta_matrix_similarity")
  results <- meta_collect_results(analyzers, result_index, allowed_calc_types=c("pairwise_distance", "pairwise_correlation"),
                                  function_name="meta_matrix_similarity")
  comparisons <- meta_build_comparisons(results, compare_label="matrix_similarity")
  out <- lapply(seq_along(comparisons), function(i) meta_compare_matrix_similarity(comparisons[[i]])) %>%
    setNames(names(comparisons))
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

#' Meta comparison for regression outputs
#' @description Compare regression outputs within or across results.
#' @param object a S4 object of class metime_analyser or a list of two metime_analyser objects
#' @param result_index character/numeric input for results. If NULL, all matching results are used.
#' @param name a character input to set the name of the results
#' @return An S4 object of class meta_analyser with the compared results and meta results
#' @export
setGeneric("meta_regression", function(object, result_index=NULL, name="meta_regression_1") standardGeneric("meta_regression"))
meta_regression_impl <- function(object, result_index=NULL, name="meta_regression_1") {
  analyzers <- meta_unpack_analyzers(object, function_name="meta_regression")
  results <- meta_collect_results(analyzers, result_index, allowed_calc_types="regression", function_name="meta_regression")
  comparisons <- meta_build_regression_comparisons(results)
  out <- lapply(seq_along(comparisons), function(i) meta_compare_regression(comparisons[[i]])) %>%
    setNames(names(comparisons))
  return(meta_make_analyser(analyzers, results, out, calc_type="meta_regression",
                            calc_info=names(out), name=name, function_name="meta_regression",
                            params=list(result_index=result_index)))
}

setMethod("meta_regression", "metime_analyser", function(object, result_index=NULL, name="meta_regression_1") {
  meta_regression_impl(object, result_index, name)
})

setMethod("meta_regression", "list", function(object, result_index=NULL, name="meta_regression_1") {
  meta_regression_impl(object, result_index, name)
})

#' Meta comparison for network overlap
#' @description Compare network edge overlap within or across results.
#' @param object a S4 object of class metime_analyser or a list of two metime_analyser objects
#' @param result_index character/numeric input for results. If NULL, all matching results are used.
#' @param name a character input to set the name of the results
#' @return An S4 object of class meta_analyser with the compared results and meta results
#' @export
setGeneric("meta_network_overlap", function(object, result_index=NULL, name="meta_network_overlap_1") standardGeneric("meta_network_overlap"))
meta_network_overlap_impl <- function(object, result_index=NULL, name="meta_network_overlap_1") {
  analyzers <- meta_unpack_analyzers(object, function_name="meta_network_overlap")
  results <- meta_collect_results(analyzers, result_index, allowed_calc_types=c("genenet_ggm", "multibipartite_ggm", "temporal_network"),
                                  function_name="meta_network_overlap")
  comparisons <- meta_build_comparisons(results, compare_label="network_overlap", allow_network_mismatch=TRUE)
  out <- lapply(seq_along(comparisons), function(i) meta_compare_network(comparisons[[i]])) %>%
    setNames(names(comparisons))
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

#' Meta comparison for feature overlap
#' @description Compare feature selection outputs within or across results.
#' @param object a S4 object of class metime_analyser or a list of two metime_analyser objects
#' @param result_index character/numeric input for results. If NULL, all matching results are used.
#' @param name a character input to set the name of the results
#' @return An S4 object of class meta_analyser with the compared results and meta results
#' @export
setGeneric("meta_feature_overlap", function(object, result_index=NULL, name="meta_feature_overlap_1") standardGeneric("meta_feature_overlap"))
meta_feature_overlap_impl <- function(object, result_index=NULL, name="meta_feature_overlap_1") {
  analyzers <- meta_unpack_analyzers(object, function_name="meta_feature_overlap")
  results <- meta_collect_results(analyzers, result_index, allowed_calc_types="feature_selection", function_name="meta_feature_overlap")
  comparisons <- meta_build_comparisons(results, compare_label="feature_overlap")
  out <- lapply(seq_along(comparisons), function(i) meta_compare_feature_overlap(comparisons[[i]])) %>%
    setNames(names(comparisons))
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

meta_unpack_analyzers <- function(object, function_name) {
  if (inherits(object, "metime_analyser")) {
    return(list(object))
  }
  if (is.list(object) && length(object) == 2 && all(vapply(object, inherits, logical(1), "metime_analyser"))) {
    return(object)
  }
  stop(paste0(function_name, "(): object must be a metime_analyser or a list of two metime_analyser objects."))
}

meta_collect_results <- function(analyzers, result_index, allowed_calc_types, function_name) {
  results <- lapply(seq_along(analyzers), function(i) {
    this_index <- result_index
    if (is.list(result_index) && length(result_index) == length(analyzers)) {
      this_index <- result_index[[i]]
    }
    meta_resolve_results(analyzers[[i]], this_index, allowed_calc_types, function_name)
  })
  results <- lapply(results, meta_normalize_result_names)
  return(results)
}

meta_resolve_results <- function(analyzer, result_index, allowed_calc_types, function_name) {
  results <- analyzer@results
  if (length(results) == 0) {
    stop(paste0(function_name, "(): no results available"))
  }
  if (is.null(result_index)) {
    result_index <- seq_along(results)
  } else if (is.character(result_index)) {
    result_index <- match(result_index, names(results))
  }
  if (any(is.na(result_index))) {
    stop(paste0(function_name, "(): result_index not found"))
  }
  selected <- results[result_index]
  if (!is.null(allowed_calc_types)) {
    invalid <- vapply(selected, function(res) {
      !all(res$information$calc_type %in% allowed_calc_types)
    }, logical(1))
    if (any(invalid)) {
      stop(paste0(function_name, "(): result calc_type does not match required types"))
    }
  }
  return(selected)
}

meta_extract_stratifications <- function(result) {
  if (is.null(result$functions_applied) || length(result$functions_applied) == 0) {
    return(NA_character_)
  }
  last_fun <- result$functions_applied[length(result$functions_applied)]
  match <- regmatches(last_fun, regexpr("stratifications\\s*=\\s*[^,\\)]*", last_fun))
  if (length(match) == 0) {
    return(NA_character_)
  }
  gsub("stratifications\\s*=\\s*", "", match)
}

meta_warn_stratification_mismatch <- function(result1, result2) {
  strat1 <- meta_extract_stratifications(result1)
  strat2 <- meta_extract_stratifications(result2)
  if (!is.na(strat1) && !is.na(strat2) && !identical(strat1, strat2)) {
    warning(paste0("Stratifications differ between results: ", strat1, " vs ", strat2))
  }
}

meta_build_comparisons <- function(results, compare_label, allow_network_mismatch=FALSE) {
  if (length(results) == 1) {
    return(meta_build_comparisons_within(results[[1]], compare_label))
  }
  if (length(results) == 2) {
    return(meta_build_comparisons_across(results[[1]], results[[2]], allow_network_mismatch))
  }
  stop(paste0("Comparison for ", compare_label, " requires one analyser or two analysers."))
}

meta_build_comparisons_within <- function(result, compare_label) {
  result <- meta_normalize_plot_names(result)
  if (length(result$plot_data) <= 1) {
    stop(paste0("Within result comparison is not possible for ", compare_label, ": only one plot_data entry."))
  }
  plot_names <- names(result$plot_data)
  if (is.null(plot_names) || any(plot_names == "")) {
    plot_names <- paste0("plot_", seq_along(result$plot_data))
  }
  plot_indices <- seq_along(result$plot_data)
  combn_idx <- combn(plot_indices, 2, simplify=FALSE)
  comparisons <- lapply(combn_idx, function(pair) {
    list(
      result1=result,
      result2=result,
      label1=plot_names[pair[1]],
      label2=plot_names[pair[2]],
      plot1=result$plot_data[[pair[1]]],
      plot2=result$plot_data[[pair[2]]]
    )
  })
  names(comparisons) <- vapply(combn_idx, function(pair) paste(plot_names[pair], collapse="__"), character(1))
  comparisons
}

meta_build_comparisons_across <- function(results1, results2, allow_network_mismatch=FALSE) {
  shared_results <- intersect(names(results1), names(results2))
  if (length(shared_results) == 0) {
    stop("These results are not comparable: no shared result names.")
  }
  comparisons <- list()
  for (res_name in shared_results) {
    res1 <- meta_normalize_plot_names(results1[[res_name]])
    res2 <- meta_normalize_plot_names(results2[[res_name]])
    if (length(unique(res1$information$calc_type)) > 1 || length(unique(res2$information$calc_type)) > 1) {
      stop("Across result comparison is not feasible because calc_type contains multiple values.")
    }
    if (!allow_network_mismatch && length(unique(c(res1$information$calc_type, res2$information$calc_type))) > 1) {
      stop("Across result comparison is not feasible because calc_type differs.")
    }
    if (allow_network_mismatch) {
      mismatch <- length(unique(c(res1$information$calc_type, res2$information$calc_type))) > 1
      if (mismatch && all(c("genenet_ggm", "multibipartite_ggm") %in% unique(c(res1$information$calc_type, res2$information$calc_type)))) {
        warning("Comparing genenet_ggm to multibipartite_ggm results; interpret overlaps with caution.")
      } else if (mismatch) {
        stop("Across result comparison is not feasible because calc_type differs.")
      }
    }
    meta_warn_stratification_mismatch(res1, res2)
    shared <- intersect(names(res1$plot_data), names(res2$plot_data))
    if (length(shared) == 0) {
      stop("These results are not comparable: no shared plot_data names.")
    }
    for (name in shared) {
      comparisons[[paste(res_name, name, sep="__")]] <- list(
        result1=res1,
        result2=res2,
        label1=res_name,
        label2=res_name,
        plot1=res1$plot_data[[name]],
        plot2=res2$plot_data[[name]]
      )
    }
  }
  comparisons
}

meta_build_regression_comparisons <- function(results) {
  if (length(results) == 1) {
    return(meta_build_comparisons_within(results[[1]], "regression"))
  }
  if (length(results) == 2) {
    return(meta_build_regression_comparisons_across(results[[1]], results[[2]]))
  }
  stop("Regression comparison requires one analyser or two analysers.")
}

meta_build_regression_comparisons_across <- function(results1, results2) {
  shared_results <- intersect(names(results1), names(results2))
  if (length(shared_results) == 0) {
    stop("These regression results are not comparable: no shared result names.")
  }
  comparisons <- list()
  for (res_name in shared_results) {
    res1 <- results1[[res_name]]
    res2 <- results2[[res_name]]
    meta_warn_stratification_mismatch(res1, res2)
    shared <- intersect(res1$information$calc_info, res2$information$calc_info)
    if (length(shared) == 0) {
      stop("These regression results are not comparable: no shared calc_info entries.")
    }
    for (info in shared) {
      idx1 <- which(res1$information$calc_info == info)
      idx2 <- which(res2$information$calc_info == info)
      if (length(idx1) == 0 || length(idx2) == 0) {
        next
      }
      comparisons[[paste(res_name, info, sep="__")]] <- list(
        result1=res1,
        result2=res2,
        label1=res_name,
        label2=res_name,
        plot1=res1$plot_data[[idx1[1]]],
        plot2=res2$plot_data[[idx2[1]]]
      )
    }
  }
  comparisons
}

meta_make_analyser <- function(analyzers, results, out, calc_type, calc_info, name, function_name, params) {
  base <- analyzers[[1]]
  source_results <- meta_merge_source_results(results)
  meta_results <- list()
  meta_results[[name]] <- list(
    functions_applied=list(meta_format_function_info(function_name, params)),
    plot_data=out,
    information=list(calc_type=rep(calc_type, each=length(out)), calc_info=calc_info),
    plots=list()
  )
  meta_object <- new("meta_analyser",
                     list_of_data=base@list_of_data,
                     list_of_col_data=base@list_of_col_data,
                     list_of_row_data=base@list_of_row_data,
                     annotations=base@annotations,
                     results=source_results,
                     meta_results=meta_results)
  return(meta_object)
}

meta_merge_source_results <- function(results) {
  if (length(results) == 1) {
    return(results[[1]])
  }
  names1 <- names(results[[1]])
  names2 <- names(results[[2]])
  if (any(names1 %in% names2)) {
    names(results[[1]]) <- paste0("analyzer1_", names1)
    names(results[[2]]) <- paste0("analyzer2_", names2)
  }
  c(results[[1]], results[[2]])
}

meta_format_function_info <- function(function_name, params) {
  param_str <- paste0(names(params), "=", params, collapse=", ")
  paste0(function_name, "(", param_str, ")")
}

meta_build_conservation_comparisons <- function(results) {
  meta_build_comparisons(results, compare_label="conservation")
}

meta_normalize_result_names <- function(results) {
  if (is.null(names(results)) || any(names(results) == "")) {
    names(results) <- paste0("result_", seq_along(results))
  }
  results
}

meta_normalize_plot_names <- function(result) {
  if (is.null(names(result$plot_data)) || any(names(result$plot_data) == "")) {
    names(result$plot_data) <- paste0("plot_", seq_along(result$plot_data))
  }
  result
}

meta_compare_conservation <- function(comp) {
  df1 <- comp$plot1
  df2 <- comp$plot2
  key1 <- if ("id" %in% names(df1)) df1$id else rownames(df1)
  key2 <- if ("id" %in% names(df2)) df2$id else rownames(df2)
  common <- intersect(key1, key2)
  if (length(common) == 0) {
    stop("These results are not comparable: no overlapping metabolites for conservation.")
  }
  df1 <- df1[key1 %in% common, , drop=FALSE]
  df2 <- df2[key2 %in% common, , drop=FALSE]
  df1 <- df1[match(common, key1), , drop=FALSE]
  df2 <- df2[match(common, key2), , drop=FALSE]
  if (!all(c("ci") %in% names(df1)) || !all(c("ci") %in% names(df2))) {
    stop("Conservation comparison requires a 'ci' column.")
  }
  cor_val <- stats::cor(df1$ci, df2$ci, use="complete.obs")
  data.frame(
    result1=comp$label1,
    result2=comp$label2,
    n_common=length(common),
    ci_correlation=cor_val,
    stringsAsFactors=FALSE
  )
}

meta_compare_matrix_similarity <- function(comp) {
  df1 <- comp$plot1
  df2 <- comp$plot2
  required <- c("row", "column", "dist")
  if (!all(required %in% names(df1)) || !all(required %in% names(df2))) {
    stop("Matrix similarity comparison requires row/column/dist columns.")
  }
  df1$key <- paste(df1$row, df1$column, sep="__")
  df2$key <- paste(df2$row, df2$column, sep="__")
  common <- intersect(df1$key, df2$key)
  if (length(common) == 0) {
    stop("These results are not comparable: no overlapping matrix pairs.")
  }
  df1 <- df1[df1$key %in% common, , drop=FALSE]
  df2 <- df2[df2$key %in% common, , drop=FALSE]
  df1 <- df1[match(common, df1$key), , drop=FALSE]
  df2 <- df2[match(common, df2$key), , drop=FALSE]
  cor_val <- stats::cor(df1$dist, df2$dist, use="complete.obs")
  data.frame(
    result1=comp$label1,
    result2=comp$label2,
    n_common=length(common),
    dist_correlation=cor_val,
    stringsAsFactors=FALSE
  )
}

meta_compare_regression <- function(comp) {
  df1 <- comp$plot1
  df2 <- comp$plot2
  keys <- intersect(names(df1), names(df2))
  join_cols <- intersect(keys, c("met", "trait"))
  if (length(join_cols) == 0) {
    df1$key <- rownames(df1)
    df2$key <- rownames(df2)
  } else {
    df1$key <- apply(df1[, join_cols, drop=FALSE], 1, paste, collapse="__")
    df2$key <- apply(df2[, join_cols, drop=FALSE], 1, paste, collapse="__")
  }
  common <- intersect(df1$key, df2$key)
  if (length(common) == 0) {
    stop("These regression results are not comparable: no overlapping rows.")
  }
  df1 <- df1[df1$key %in% common, , drop=FALSE]
  df2 <- df2[df2$key %in% common, , drop=FALSE]
  df1 <- df1[match(common, df1$key), , drop=FALSE]
  df2 <- df2[match(common, df2$key), , drop=FALSE]
  if (!all(c("beta") %in% names(df1)) || !all(c("beta") %in% names(df2))) {
    stop("Regression comparison requires a 'beta' column.")
  }
  cor_beta <- stats::cor(df1$beta, df2$beta, use="complete.obs")
  sign_agreement <- mean(sign(df1$beta) == sign(df2$beta), na.rm=TRUE)
  data.frame(
    result1=comp$label1,
    result2=comp$label2,
    n_common=length(common),
    beta_correlation=cor_beta,
    sign_agreement=sign_agreement,
    stringsAsFactors=FALSE
  )
}

meta_compare_network <- function(comp) {
  edge1 <- meta_extract_edges(comp$plot1)
  edge2 <- meta_extract_edges(comp$plot2)
  if (!all(c("node1", "node2") %in% names(edge1)) || !all(c("node1", "node2") %in% names(edge2))) {
    stop("Network overlap comparison requires node1/node2 columns.")
  }
  type1 <- unique(comp$result1$information$calc_type)
  type2 <- unique(comp$result2$information$calc_type)
  directed <- length(intersect(c(type1, type2), "temporal_network")) >= 1
  edge_id1 <- meta_make_edge_ids(edge1, directed=directed)
  edge_id2 <- meta_make_edge_ids(edge2, directed=directed)
  common <- intersect(edge_id1, edge_id2)
  if (length(common) == 0) {
    stop("These results are not comparable: no overlapping edges.")
  }
  union_edges <- length(unique(c(edge_id1, edge_id2)))
  data.frame(
    result1=comp$label1,
    result2=comp$label2,
    n_common=length(common),
    jaccard=length(common) / union_edges,
    stringsAsFactors=FALSE
  )
}

meta_compare_feature_overlap <- function(comp) {
  df1 <- comp$plot1
  df2 <- comp$plot2
  if (!all(c("met", "selected") %in% names(df1)) || !all(c("met", "selected") %in% names(df2))) {
    stop("Feature overlap comparison requires met and selected columns.")
  }
  selected1 <- df1$met[df1$selected %in% TRUE]
  selected2 <- df2$met[df2$selected %in% TRUE]
  common <- intersect(selected1, selected2)
  if (length(common) == 0) {
    stop("These results are not comparable: no overlapping selected features.")
  }
  union_features <- length(unique(c(selected1, selected2)))
  data.frame(
    result1=comp$label1,
    result2=comp$label2,
    n_common=length(common),
    jaccard=length(common) / union_features,
    stringsAsFactors=FALSE
  )
}

meta_extract_edges <- function(plot_entry) {
  if (is.list(plot_entry) && !is.data.frame(plot_entry) && !is.null(plot_entry$edge)) {
    return(plot_entry$edge)
  }
  if (is.data.frame(plot_entry)) {
    return(plot_entry)
  }
  stop("Network overlap comparison requires an edge list.")
}

meta_make_edge_ids <- function(edge_df, directed=FALSE) {
  if (directed) {
    return(paste(edge_df$node1, edge_df$node2, sep="__"))
  }
  node_a <- pmin(edge_df$node1, edge_df$node2)
  node_b <- pmax(edge_df$node1, edge_df$node2)
  paste(node_a, node_b, sep="__")
}
