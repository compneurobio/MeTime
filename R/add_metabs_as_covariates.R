#' Function to add metabolites as covariates for network construction
#' @description Method applied on metime_analyser object to add other metabolite data to a certain dataset
#' @param object A S4 object of class metime_analyser
#' @param which_data Dataset to which the metab data is to be added(please note that this a single character)
#' @param which_metabs list of names of metabs and name of the list represents the dataset from which 
#' the metabs are to be acquired. eg: which_metabs=list(nmr_data=c("metab1", "metab2"), lipid_data=c(""))
#' @return S4 object with metabs added for GGM to another dataset
#' @export
setGeneric("add_metabs_as_covariates", function(object, which_data, which_metabs) standardGeneric("add_metabs_as_covariates"))
setMethod("add_metabs_as_covariates", "metime_analyser", function(object, which_data, which_metabs) {
			data <- object@list_of_data[[which_data]]
			list_of_metabs <- list()
			list_of_metabs <- lapply(seq_along(which_metabs), function(i) {
				dat <- object@list_of_data[[names(which_metabs)[i]]]
				dat <- dat[ ,which_metabs[[i]]]
				dat <- dat[order(rownames(dat)), ]
				return(dat)
			})
			if(length(which_data) > 1) {
				metab_matrix <- do.call(cbind, list_of_metabs)
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


