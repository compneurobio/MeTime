#' Validate a metime_analyser object
#' @description Check consistency of list_of_data, row_data, and col_data.
#' @param object a S4 object of class "metime_analyser".
#' @param which_data character vector of dataset names to validate. If NULL, validates all datasets.
#' @return A list with dataset-level issues and a summary table.
#' @export
setGeneric("validate_metime_analyser", function(object, which_data=NULL) standardGeneric("validate_metime_analyser"))
setMethod("validate_metime_analyser", "metime_analyser", function(object, which_data=NULL) {
  data_names <- names(object@list_of_data)
  if(is.null(which_data)) which_data <- data_names
  missing <- setdiff(which_data, data_names)
  if(length(missing) > 0) warning("Datasets not found: ", paste(missing, collapse=", "))
  which_data <- intersect(which_data, data_names)

  issues <- lapply(which_data, function(dataset) {
    data <- object@list_of_data[[dataset]]
    row_data <- object@list_of_row_data[[dataset]]
    col_data <- object@list_of_col_data[[dataset]]
    dataset_issues <- c()

    if(is.null(data)) dataset_issues <- c(dataset_issues, "data missing")
    if(is.null(row_data)) dataset_issues <- c(dataset_issues, "row_data missing")
    if(is.null(col_data)) dataset_issues <- c(dataset_issues, "col_data missing")

    if(!is.null(data) && !is.null(row_data)) {
      if(!all(rownames(data) %in% row_data$id)) dataset_issues <- c(dataset_issues, "row_data ids do not match rownames")
      if(nrow(data) != nrow(row_data)) dataset_issues <- c(dataset_issues, "row_data rows do not match data rows")
      if(!all(c("id","subject","time") %in% names(row_data))) dataset_issues <- c(dataset_issues, "row_data missing required columns")
    }

    if(!is.null(data) && !is.null(col_data)) {
      if("id" %in% names(col_data)) {
        if(!all(colnames(data) %in% col_data$id)) dataset_issues <- c(dataset_issues, "col_data ids do not match colnames")
      }
      if(ncol(data) != nrow(col_data)) dataset_issues <- c(dataset_issues, "col_data rows do not match data columns")
    }

    dataset_issues
  })
  names(issues) <- which_data

  summary <- dplyr::tibble(
    dataset = which_data,
    issue_count = vapply(issues, length, integer(1)),
    ok = vapply(issues, function(x) length(x) == 0, logical(1))
  )

  list(ok = all(summary$ok), issues = issues, summary = summary)
})

#' Summarize datasets in a metime_analyser object
#' @description Returns dataset size and missingness summaries.
#' @param object a S4 object of class "metime_analyser".
#' @param which_data character vector of dataset names to summarize. If NULL, summarizes all datasets.
#' @return A data.frame with summary statistics for each dataset.
#' @export
setGeneric("summarize_dataset", function(object, which_data=NULL) standardGeneric("summarize_dataset"))
setMethod("summarize_dataset", "metime_analyser", function(object, which_data=NULL) {
  data_names <- names(object@list_of_data)
  if(is.null(which_data)) which_data <- data_names
  which_data <- intersect(which_data, data_names)

  summaries <- lapply(which_data, function(dataset) {
    data <- object@list_of_data[[dataset]]
    row_data <- object@list_of_row_data[[dataset]]

    total_values <- length(data)
    missing_values <- sum(is.na(data))
    timepoints <- if(!is.null(row_data) && "time" %in% names(row_data)) length(unique(row_data$time)) else NA_integer_
    subjects <- if(!is.null(row_data) && "subject" %in% names(row_data)) length(unique(row_data$subject)) else NA_integer_

    dplyr::tibble(
      dataset = dataset,
      n_samples = nrow(data),
      n_features = ncol(data),
      missing_values = missing_values,
      missing_pct = ifelse(total_values > 0, missing_values / total_values, NA_real_),
      n_timepoints = timepoints,
      n_subjects = subjects
    )
  })

  dplyr::bind_rows(summaries)
})

#' Plot missingness heatmap for a dataset
#' @description Creates a heatmap of missing values for a dataset.
#' @param object a S4 object of class "metime_analyser".
#' @param which_data a single dataset name.
#' @param max_features maximum number of features to include in the plot.
#' @return A ggplot object.
#' @export
setGeneric("plot_missingness", function(object, which_data, max_features=100) standardGeneric("plot_missingness"))
setMethod("plot_missingness", "metime_analyser", function(object, which_data, max_features=100) {
  if(length(which_data) != 1) stop("plot_missingness expects a single dataset name")
  if(!which_data %in% names(object@list_of_data)) stop("Dataset not found in object")

  data <- object@list_of_data[[which_data]]
  if(ncol(data) > max_features) {
    data <- data[, seq_len(max_features), drop=FALSE]
  }

  missing <- as.data.frame(is.na(data))
  missing$sample <- rownames(missing)
  missing_long <- tidyr::pivot_longer(missing, cols = -sample, names_to = "feature", values_to = "missing")

  ggplot2::ggplot(missing_long, ggplot2::aes(x = feature, y = sample, fill = missing)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#2c7fb8", "TRUE" = "#d95f0e")) +
    ggplot2::labs(title = paste("Missingness:", which_data), x = "Feature", y = "Sample") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
})

#' Filter timepoints from datasets
#' @description Removes samples that are not in the specified timepoints.
#' @param object a S4 object of class "metime_analyser".
#' @param timepoints character or numeric vector of timepoints to keep.
#' @param which_data character vector of dataset names to filter. If NULL, filters all datasets.
#' @return A modified metime_analyser object.
#' @export
setGeneric("mod_filter_timepoints", function(object, timepoints, which_data=NULL) standardGeneric("mod_filter_timepoints"))
setMethod("mod_filter_timepoints", "metime_analyser", function(object, timepoints, which_data=NULL) {
  data_names <- names(object@list_of_data)
  if(is.null(which_data)) which_data <- data_names
  which_data <- intersect(which_data, data_names)

  for(dataset in which_data) {
    row_data <- object@list_of_row_data[[dataset]]
    if(is.null(row_data) || !"time" %in% names(row_data)) {
      warning("time column not found in row_data for dataset: ", dataset)
      next
    }
    keep <- row_data$time %in% timepoints
    keep_ids <- rownames(row_data)[keep]
    object@list_of_row_data[[dataset]] <- row_data[keep, , drop=FALSE]
    object@list_of_data[[dataset]] <- object@list_of_data[[dataset]][keep_ids, , drop=FALSE]
  }

  out <- add_function_info(object=object, function_name="mod_filter_timepoints", params=list(timepoints=timepoints, which_data=which_data))
  return(out)
})

#' Filter subjects from datasets
#' @description Removes samples that are not in the specified subjects.
#' @param object a S4 object of class "metime_analyser".
#' @param subjects character vector of subjects to keep.
#' @param which_data character vector of dataset names to filter. If NULL, filters all datasets.
#' @return A modified metime_analyser object.
#' @export
setGeneric("mod_filter_subjects", function(object, subjects, which_data=NULL) standardGeneric("mod_filter_subjects"))
setMethod("mod_filter_subjects", "metime_analyser", function(object, subjects, which_data=NULL) {
  data_names <- names(object@list_of_data)
  if(is.null(which_data)) which_data <- data_names
  which_data <- intersect(which_data, data_names)

  for(dataset in which_data) {
    row_data <- object@list_of_row_data[[dataset]]
    if(is.null(row_data) || !"subject" %in% names(row_data)) {
      warning("subject column not found in row_data for dataset: ", dataset)
      next
    }
    keep <- row_data$subject %in% subjects
    keep_ids <- rownames(row_data)[keep]
    object@list_of_row_data[[dataset]] <- row_data[keep, , drop=FALSE]
    object@list_of_data[[dataset]] <- object@list_of_data[[dataset]][keep_ids, , drop=FALSE]
  }

  out <- add_function_info(object=object, function_name="mod_filter_subjects", params=list(subjects=subjects, which_data=which_data))
  return(out)
})

#' Filter features by missingness
#' @description Removes features with missingness above a threshold.
#' @param object a S4 object of class "metime_analyser".
#' @param threshold numeric between 0 and 1 defining maximum allowed missingness.
#' @param which_data character vector of dataset names to filter. If NULL, filters all datasets.
#' @return A modified metime_analyser object.
#' @export
setGeneric("mod_filter_features_by_missingness", function(object, threshold=0.2, which_data=NULL) standardGeneric("mod_filter_features_by_missingness"))
setMethod("mod_filter_features_by_missingness", "metime_analyser", function(object, threshold=0.2, which_data=NULL) {
  data_names <- names(object@list_of_data)
  if(is.null(which_data)) which_data <- data_names
  which_data <- intersect(which_data, data_names)

  for(dataset in which_data) {
    data <- object@list_of_data[[dataset]]
    if(is.null(data)) next
    missing_rate <- colMeans(is.na(data))
    keep <- missing_rate <= threshold
    data <- data[, keep, drop=FALSE]
    object@list_of_data[[dataset]] <- data

    col_data <- object@list_of_col_data[[dataset]]
    if(!is.null(col_data)) {
      col_ids <- if("id" %in% names(col_data)) col_data$id else rownames(col_data)
      object@list_of_col_data[[dataset]] <- col_data[col_ids %in% colnames(data), , drop=FALSE]
    }
  }

  out <- add_function_info(object=object, function_name="mod_filter_features_by_missingness", params=list(threshold=threshold, which_data=which_data))
  return(out)
})

#' Filter samples by missingness
#' @description Removes samples with missingness above a threshold.
#' @param object a S4 object of class "metime_analyser".
#' @param threshold numeric between 0 and 1 defining maximum allowed missingness.
#' @param which_data character vector of dataset names to filter. If NULL, filters all datasets.
#' @return A modified metime_analyser object.
#' @export
setGeneric("mod_filter_samples_by_missingness", function(object, threshold=0.2, which_data=NULL) standardGeneric("mod_filter_samples_by_missingness"))
setMethod("mod_filter_samples_by_missingness", "metime_analyser", function(object, threshold=0.2, which_data=NULL) {
  data_names <- names(object@list_of_data)
  if(is.null(which_data)) which_data <- data_names
  which_data <- intersect(which_data, data_names)

  for(dataset in which_data) {
    data <- object@list_of_data[[dataset]]
    if(is.null(data)) next
    missing_rate <- rowMeans(is.na(data))
    keep <- missing_rate <= threshold
    data <- data[keep, , drop=FALSE]
    object@list_of_data[[dataset]] <- data

    row_data <- object@list_of_row_data[[dataset]]
    if(!is.null(row_data)) {
      object@list_of_row_data[[dataset]] <- row_data[rownames(row_data) %in% rownames(data), , drop=FALSE]
    }
  }

  out <- add_function_info(object=object, function_name="mod_filter_samples_by_missingness", params=list(threshold=threshold, which_data=which_data))
  return(out)
})

#' Filter features by variance
#' @description Removes features with variance below a threshold.
#' @param object a S4 object of class "metime_analyser".
#' @param min_variance numeric defining minimum variance required to keep a feature.
#' @param which_data character vector of dataset names to filter. If NULL, filters all datasets.
#' @return A modified metime_analyser object.
#' @export
setGeneric("mod_filter_features_by_variance", function(object, min_variance=0, which_data=NULL) standardGeneric("mod_filter_features_by_variance"))
setMethod("mod_filter_features_by_variance", "metime_analyser", function(object, min_variance=0, which_data=NULL) {
  data_names <- names(object@list_of_data)
  if(is.null(which_data)) which_data <- data_names
  which_data <- intersect(which_data, data_names)

  for(dataset in which_data) {
    data <- object@list_of_data[[dataset]]
    if(is.null(data)) next
    variances <- apply(data, 2, stats::var, na.rm=TRUE)
    keep <- variances >= min_variance
    data <- data[, keep, drop=FALSE]
    object@list_of_data[[dataset]] <- data

    col_data <- object@list_of_col_data[[dataset]]
    if(!is.null(col_data)) {
      col_ids <- if("id" %in% names(col_data)) col_data$id else rownames(col_data)
      object@list_of_col_data[[dataset]] <- col_data[col_ids %in% colnames(data), , drop=FALSE]
    }
  }

  out <- add_function_info(object=object, function_name="mod_filter_features_by_variance", params=list(min_variance=min_variance, which_data=which_data))
  return(out)
})

#' Calculate baseline change per subject
#' @description Computes baseline-adjusted values for each sample.
#' @param object a S4 object of class "metime_analyser".
#' @param which_data dataset name to use.
#' @param baseline_timepoint timepoint used as baseline.
#' @param method "diff" for subtraction or "ratio" for division.
#' @param name name for the results entry.
#' @return A metime_analyser object with results appended.
#' @export
setGeneric("calc_baseline_change", function(object, which_data, baseline_timepoint, method="diff", name=NULL) standardGeneric("calc_baseline_change"))
setMethod("calc_baseline_change", "metime_analyser", function(object, which_data, baseline_timepoint, method="diff", name=NULL) {
  if(!which_data %in% names(object@list_of_data)) stop("Dataset not found in object")
  data <- object@list_of_data[[which_data]]
  row_data <- object@list_of_row_data[[which_data]]

  if(is.null(row_data) || !all(c("subject", "time") %in% names(row_data))) {
    stop("row_data must contain subject and time columns")
  }

  baseline_rows <- row_data$time %in% baseline_timepoint
  if(!any(baseline_rows)) {
    warning("No baseline rows found for timepoint: ", baseline_timepoint)
    return(object)
  }

  baseline_data <- data[baseline_rows, , drop=FALSE]
  baseline_subjects <- row_data$subject[baseline_rows]
  baseline_index <- match(row_data$subject, baseline_subjects)
  baseline_matrix <- baseline_data[baseline_index, , drop=FALSE]

  if(method == "diff") {
    delta <- data - baseline_matrix
  } else if(method == "ratio") {
    delta <- data / baseline_matrix
  } else {
    stop("method must be 'diff' or 'ratio'")
  }

  delta <- as.data.frame(delta)
  rownames(delta) <- rownames(data)
  result_name <- if(is.null(name)) paste0(which_data, "_baseline_change") else name

  out <- get_make_results(object=object, data=list(delta), metadata=NULL,
    calc_type="baseline_change",
    calc_info=paste("baseline_timepoint:", baseline_timepoint, "method:", method),
    name=result_name)
  out <- add_function_info(object=out, function_name="calc_baseline_change",
    params=list(which_data=which_data, baseline_timepoint=baseline_timepoint, method=method))
  return(out)
})

#' Calculate feature time trends
#' @description Fits a linear trend per feature over time.
#' @param object a S4 object of class "metime_analyser".
#' @param which_data dataset name to use.
#' @param name name for the results entry.
#' @return A metime_analyser object with results appended.
#' @export
setGeneric("calc_time_trend", function(object, which_data, name=NULL) standardGeneric("calc_time_trend"))
setMethod("calc_time_trend", "metime_analyser", function(object, which_data, name=NULL) {
  if(!which_data %in% names(object@list_of_data)) stop("Dataset not found in object")
  data <- object@list_of_data[[which_data]]
  row_data <- object@list_of_row_data[[which_data]]

  if(is.null(row_data) || !"time" %in% names(row_data)) {
    stop("row_data must contain a time column")
  }

  time_num <- suppressWarnings(as.numeric(row_data$time))
  if(all(is.na(time_num))) stop("time column could not be coerced to numeric")

  results <- lapply(seq_len(ncol(data)), function(i) {
    feature <- colnames(data)[i]
    y <- data[, i]
    valid <- is.finite(time_num) & !is.na(y)
    if(sum(valid) < 3) {
      return(dplyr::tibble(feature = feature, slope = NA_real_, pval = NA_real_, n = sum(valid)))
    }
    model_data <- data.frame(time = time_num[valid], y = y[valid])
    model <- stats::lm(y ~ time, data = model_data)
    slope <- stats::coef(model)[2]
    pval <- summary(model)$coefficients["time", "Pr(>|t|)"]
    dplyr::tibble(feature = feature, slope = slope, pval = pval, n = sum(valid))
  })

  results <- dplyr::bind_rows(results)
  result_name <- if(is.null(name)) paste0(which_data, "_time_trend") else name

  out <- get_make_results(object=object, data=list(results), metadata=NULL,
    calc_type="time_trend",
    calc_info="Linear trend per feature over time",
    name=result_name)
  out <- add_function_info(object=out, function_name="calc_time_trend",
    params=list(which_data=which_data))
  return(out)
})
