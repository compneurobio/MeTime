
#' Modification (mod) function to merge data and row_data (sample info) partially or completely. 
#' @description a modification function that merges partial or complete data and row_data of a dataset.
#' If you want to use certain columns for the rest of the analysis properly for annotating plots then use mod_mutate 
#' function to rename the colnames and rownames (if applicable) to maintain consistency in the merged dataset 
#' @param object an S4 object of class metime_analyser.
#' @param which_data a character defining which data and row_data should be merged. Has to contain only one value.
#' @param cols_list A list of named character vectors to merge dataset. 
#' example: which_data = "dataset1"; cols_list = list(data=c(...), row_data=c(...)) or list(c(...), c(...))
#' If you want full data set cols_list=list(data=NULL, row_data=NULL) which is the default setting
#' @param name a character to define the name of a new dataset.
#' @return a new S4 object of class metime_analyser with the  new merged dataset appended to it
#' @seealso [mod_merge_results], [mod_merge_data]
#' @export 
#' 
setGeneric("mod_merge_row_data_and_data", function(object, which_data, cols_list=list(data = NULL, row_data = NULL), append=FALSE, name="merged_data") standardGeneric("mod_merge_row_data_and_data")) 
setMethod("mod_merge_row_data_and_data", "metime_analyser", function(object, which_data,
                                    cols_list = list(data = NULL, row_data = NULL),
                                    append = FALSE,
                                    name = NULL) {
  if (length(which_data) != 1) stop("which_data must be a single dataset name")
  if (!which_data %in% names(object@list_of_data)) stop("Dataset not found in object")
  data <- object@list_of_data[[which_data]]
  row_data <- object@list_of_row_data[[which_data]]
  if (is.null(data) || is.null(row_data)) stop("data and row_data must be present")

  common_ids <- intersect(rownames(data), rownames(row_data))
  if (length(common_ids) == 0) stop("No matching rownames between data and row_data")
  data <- data[common_ids, , drop = FALSE]
  row_data <- row_data[common_ids, , drop = FALSE]

  if (!is.null(cols_list$data)) {
    data <- data[, intersect(cols_list$data, colnames(data)), drop = FALSE]
  }
  if (!is.null(cols_list$row_data)) {
    row_data <- row_data[, intersect(cols_list$row_data, colnames(row_data)), drop = FALSE]
  }

  drop_cols <- intersect(c("id", "subject", "time"), colnames(row_data))
  if (length(drop_cols) > 0) {
    row_data <- row_data[, setdiff(colnames(row_data), drop_cols), drop = FALSE]
  }

  merged <- cbind(data, row_data)

  if (append) {
    if (is.null(name) || nchar(name) == 0) {
      name <- paste0(which_data, "_row_merged")
    }
    if (name %in% names(object@list_of_data)) {
      stop("name already exists in object@list_of_data")
    }
    # append data to analyzer object
  
    object <- add_dataset(object=object, data=merged, col_data=object@list_of_col_data[[which_data]], row_data=object@list_of_row_data[[which_data]], name=name)
  } else {
    object@list_of_data[[which_data]] <- merged
  }
  object <- add_function_info(object=object, function_name="mod_merge_row_data_and_data", params=list(which_data=which_data,
      mutations=cols_list, name=name))
  return(object)
})