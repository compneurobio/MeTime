#' Add columns from data x to data y within a metime_analyzer.
#' @description Merge a specified columns from datasets within the S4 object to one data.
#' @param object a S4 object of the class "metime_analyzer".
#' @param which_data a character defining the dataset to which the metab data is to be added. 
#' @param which_metabs a named list of character vectors. Names of the list elements correspond to the name of dataset, the character vectors define columns to be used. eg: which_metabs=list(nmr_data=c("metab1", "metab2"), lipid_data=c("metab3", "metab4"))
#' @return object of class metime_analyser with columns from one dataset added for analysis to another dataset
#' @export
setGeneric("add_metabs_as_covariates", function(object, which_data, which_metabs) standardGeneric("add_metabs_as_covariates"))
setMethod("add_metabs_as_covariates", "metime_analyser", function(object, which_data, which_metabs) {
			data <- object@list_of_data[[which_data]]
			list_of_metabs <- lapply(seq_along(which_metabs), function(i) {
				dat <- object@list_of_data[[names(which_metabs)[i]]]
				dat <- dat[ ,which_metabs[[i]]]
				dat <- dat[order(rownames(dat)), ]
				return(dat)
			})
			if(length(which_metabs) > 1) {
				samples <- lapply(seq_along(list_of_metabs), function(a) {
						return(rownames(list_of_metabs))
					})
				common_samples <- Reduce(intersect, samples)
				metab_matrix <- lapply(seq_along(list_of_metabs), function(a) {
						return(list_of_metabs[[a]][rownames(list_of_metabs[[a]]) %in% common_samples, ])
					}) %>% do.call(what=cbind.data.frame)
			} else {
				metab_matrix <- list_of_metabs[[1]]
			}
			list_of_samples <- list(rownames(data), rownames(metab_matrix))
			common_samples <- Reduce(intersect, list_of_samples)
			data <- data[rownames(data) %in% common_samples, ]
			data <- data[order(rownames(data)),]
			metab_matrix <- metab_matrix[rownames(metab_matrix) %in% common_samples, ]
			metab_matrix <- metab_matrix[order(rownames(metab_matrix)), ]
			object@list_of_data[[which_data]] <- as.data.frame(cbind(data, metab_matrix))
			out <- object %>% add_function_info(function_name="add_metabs_as_covariates", 
				params=list(which_data=which_data, which_metabs=which_metabs))
			return(out)
	})


