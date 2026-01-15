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
    return(meta_build_comparisons_single_analyzer(results[[1]], compare_label, allow_network_mismatch))
  }
  if (length(results) == 2) {
    return(meta_build_comparisons_across(results[[1]], results[[2]], allow_network_mismatch))
  }
  stop(paste0("Comparison for ", compare_label, " requires one analyser or two analysers."))
}

meta_build_comparisons_single_analyzer <- function(results, compare_label, allow_network_mismatch=FALSE) {
  calc_groups <- meta_group_results_by_calc_type(results)
  comparisons <- list()
  for (group_name in names(calc_groups)) {
    group <- calc_groups[[group_name]]
    if (group_name %in% c("mixed", "unknown")) {
      for (res_name in names(group)) {
        comparisons <- c(comparisons,
                         meta_prefix_comparison_names(
                           meta_build_comparisons_within(group[[res_name]], compare_label),
                           res_name))
      }
      next
    }
    if (length(group) == 1) {
      res_name <- names(group)[1]
      comparisons <- c(comparisons,
                       meta_prefix_comparison_names(
                         meta_build_comparisons_within(group[[1]], compare_label),
                         res_name))
    } else {
      for (res_name in names(group)) {
        comparisons <- c(comparisons,
                         meta_prefix_comparison_names(
                           meta_build_comparisons_within(group[[res_name]], compare_label),
                           res_name))
      }
      comparisons <- c(comparisons, meta_build_comparisons_across_results(group, allow_network_mismatch))
    }
  }
  if (length(comparisons) == 0) {
    stop("No comparable results found for comparison.")
  }
  comparisons
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
  if (length(comparisons) == 0) {
    stop("These results are not comparable: no shared calc_type entries.")
  }
  comparisons
}

meta_build_regression_comparisons <- function(results) {
  if (length(results) == 1) {
    return(meta_build_regression_single_analyzer(results[[1]]))
  }
  if (length(results) == 2) {
    return(meta_build_regression_comparisons_across(results[[1]], results[[2]]))
  }
  stop("Regression comparison requires one analyser or two analysers.")
}

meta_build_regression_single_analyzer <- function(results) {
  calc_groups <- meta_group_results_by_calc_type(results)
  comparisons <- list()
  for (group_name in names(calc_groups)) {
    group <- calc_groups[[group_name]]
    if (group_name %in% c("mixed", "unknown")) {
      for (res_name in names(group)) {
        comparisons <- c(comparisons,
                         meta_prefix_comparison_names(
                           meta_build_comparisons_within(group[[res_name]], "regression"),
                           res_name))
      }
      next
    }
    if (length(group) == 1) {
      res_name <- names(group)[1]
      comparisons <- c(comparisons,
                       meta_prefix_comparison_names(
                         meta_build_comparisons_within(group[[1]], "regression"),
                         res_name))
    } else {
      for (res_name in names(group)) {
        comparisons <- c(comparisons,
                         meta_prefix_comparison_names(
                           meta_build_comparisons_within(group[[res_name]], "regression"),
                           res_name))
      }
      comparisons <- c(comparisons, meta_build_regression_comparisons_across(group, group))
    }
  }
  if (length(comparisons) == 0) {
    stop("No comparable regression results found.")
  }
  comparisons
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
  if (!isClass("meta_analyser")) {
    methods::setClass("meta_analyser",
                      slots=list(results="list", annotations="list", meta_results="list"),
                      where=globalenv())
  }
  source_results <- meta_merge_source_results(results)
  meta_results <- list()
  meta_results[[name]] <- list(
    functions_applied=list(meta_format_function_info(function_name, params)),
    plot_data=out,
    information=list(calc_type=rep(calc_type, each=length(out)), calc_info=calc_info),
    plots=list()
  )
  meta_object <- new("meta_analyser",
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
  if (is.null(result$plot_data)) {
    return(result)
  }
  if (is.null(names(result$plot_data)) || any(names(result$plot_data) == "")) {
    names(result$plot_data) <- paste0("plot_", seq_along(result$plot_data))
  }
  result
}

meta_compare_conservation <- function(comp, top_k=50) {
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
  if (!"id" %in% names(df1)) {
    df1$id <- key1
  }
  if (!"id" %in% names(df2)) {
    df2$id <- key2
  }
  out <- list()
  out[["score_cor_pearson"]] <- data.frame(
    result1=comp$label1,
    result2=comp$label2,
    n_common=length(common),
    score_correlation=stats::cor(df1$ci, df2$ci, use="complete.obs", method="pearson"),
    method="pearson",
    stringsAsFactors=FALSE
  )
  out[["score_cor_spearman"]] <- data.frame(
    result1=comp$label1,
    result2=comp$label2,
    n_common=length(common),
    score_correlation=stats::cor(df1$ci, df2$ci, use="complete.obs", method="spearman"),
    method="spearman",
    stringsAsFactors=FALSE
  )
  k <- min(top_k, length(common))
  top1 <- df1[order(df1$ci, decreasing=TRUE), , drop=FALSE]
  top2 <- df2[order(df2$ci, decreasing=TRUE), , drop=FALSE]
  top_ids1 <- top1$id[seq_len(k)]
  top_ids2 <- top2$id[seq_len(k)]
  overlap <- intersect(top_ids1, top_ids2)
  union_vals <- length(unique(c(top_ids1, top_ids2)))
  out[["topk_jaccard"]] <- data.frame(
    result1=comp$label1,
    result2=comp$label2,
    top_k=k,
    n_overlap=length(overlap),
    jaccard=ifelse(union_vals == 0, NA_real_, length(overlap) / union_vals),
    stringsAsFactors=FALSE
  )
  rank1 <- rank(-df1$ci, ties.method="average")
  rank2 <- rank(-df2$ci, ties.method="average")
  out[["rank_kendall"]] <- data.frame(
    result1=comp$label1,
    result2=comp$label2,
    n_common=length(common),
    rank_similarity=stats::cor(rank1, rank2, use="complete.obs", method="kendall"),
    method="kendall",
    stringsAsFactors=FALSE
  )
  out[["rank_spearman"]] <- data.frame(
    result1=comp$label1,
    result2=comp$label2,
    n_common=length(common),
    rank_similarity=stats::cor(rank1, rank2, use="complete.obs", method="spearman"),
    method="spearman",
    stringsAsFactors=FALSE
  )
  out
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

meta_compare_regression <- function(comp, method) {
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
  out <- list()
  if ("sign" %in% method) {
    sign_df <- data.frame(
      sign1=ifelse(df1$beta >= 0, "+", "-"),
      sign2=ifelse(df2$beta >= 0, "+", "-"),
      stringsAsFactors=FALSE
    ) %>%
      dplyr::mutate(combined=paste(sign1, sign2)) %>%
      dplyr::count(combined) %>%
      tidyr::spread(key=combined, value=n)
    out[["sign"]] <- cbind.data.frame(
      result1=comp$label1,
      result2=comp$label2,
      n_common=length(common),
      sign_df,
      stringsAsFactors=FALSE
    )
  }
  if ("het" %in% method) {
    if (!all(c("se") %in% names(df1)) || !all(c("se") %in% names(df2))) {
      stop("Heterogeneity comparison requires a 'se' column.")
    }
    diff_t <- (df1$beta - df2$beta) / sqrt(df1$se^2 + df2$se^2)
    diff_p <- 2 * stats::pnorm(-abs(diff_t))
    i_sq <- ifelse(abs(diff_t) > 1, ((abs(diff_t) - 1) / abs(diff_t)) * 100, 0)
    sig <- ifelse(diff_p <= 0.05 / length(common), "significant", "non_significant")
    het_df <- data.frame(significant=sig, stringsAsFactors=FALSE) %>%
      dplyr::count(significant) %>%
      tidyr::spread(key=significant, value=n)
    out[["het"]] <- cbind.data.frame(
      result1=comp$label1,
      result2=comp$label2,
      n_common=length(common),
      i_sq_mean=mean(i_sq, na.rm=TRUE),
      het_df,
      stringsAsFactors=FALSE
    )
  }
  if ("cor" %in% method) {
    cor_beta <- stats::cor(df1$beta, df2$beta, use="complete.obs")
    out[["cor"]] <- data.frame(
      result1=comp$label1,
      result2=comp$label2,
      n_common=length(common),
      beta_correlation=cor_beta,
      stringsAsFactors=FALSE
    )
  }
  out
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

meta_prefix_comparison_names <- function(comparisons, prefix) {
  if (length(comparisons) == 0) {
    return(comparisons)
  }
  names(comparisons) <- paste(prefix, names(comparisons), sep="__")
  comparisons
}

meta_group_results_by_calc_type <- function(results) {
  calc_types <- vapply(results, meta_result_calc_type_key, character(1))
  split(results, calc_types)
}

meta_result_calc_type_key <- function(result) {
  if (is.null(result$information$calc_type)) {
    return("unknown")
  }
  types <- unique(result$information$calc_type)
  if (length(types) == 1) {
    return(types)
  }
  "mixed"
}

meta_build_comparisons_across_results <- function(results, allow_network_mismatch=FALSE) {
  result_names <- names(results)
  if (length(result_names) < 2) {
    stop("Across result comparison requires at least two results.")
  }
  combn_pairs <- combn(result_names, 2, simplify=FALSE)
  comparisons <- list()
  for (pair in combn_pairs) {
    res1 <- meta_normalize_plot_names(results[[pair[1]]])
    res2 <- meta_normalize_plot_names(results[[pair[2]]])
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
      comparisons[[paste(pair[1], pair[2], name, sep="__")]] <- list(
        result1=res1,
        result2=res2,
        label1=pair[1],
        label2=pair[2],
        plot1=res1$plot_data[[name]],
        plot2=res2$plot_data[[name]]
      )
    }
  }
  comparisons
}