#' Function to add covariates to the dataset of interest for GGMs
#' @description adds Covariates to data matrices in metime_analyser S4 object
#' @param object object of class metime_analyser
#' @param which_data Dataset to which the covariates is to be added
#' @param covariates character vector names of covariates. 
#' @param class.ind Logical to convert factor variables into class.ind style or not
#' @param phenotype Logical. If True will extract from phenotype dataset else uses row data
#' @return S4 object with covariates added to the dataset
#' @export
setGeneric("add_phenotypes_as_covariates", function(object, which_data, covariates, class.ind=FALSE, phenotype=FALSE) standardGeneric("add_phenotypes_as_covariates"))
setMethod("add_phenotypes_as_covariates", "metime_analyser", function(object, which_data, covariates, class.ind=FALSE, phenotype=FALSE) {
			stopifnot(which_data %>% length() == 1 )
			if(phenotype) {
				phenotype_data <- object@list_of_data[[object@annotations[[1]]$phenotype]]
				phenotype_data <- phenotype_data[order(rownames(phenotype_data)),]
				list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
				covariates_list <- lapply(list_of_data, function(x) {
							x <- x[order(rownames(x)), ]
							return(phenotype_data[rownames(phenotype_data) %in% rownames(x), covariates])
					})
				for(i in 1:length(which_data)) {
					 list_of_data <- unname(list_of_data)
					 covariates_list <- unname(covariates_list)
					 list_of_data[[i]] <- as.data.frame(cbind(list_of_data[[i]], covariates_list[[i]]))
					 if(class.ind) {
					 	for(j in 1:length(covariates)) {
					 		vec <- list_of_data[[i]][ ,covariates[j]]
					 		if(class(vec) %in% "character" | class(vec) %in% "factor") {
					 			new_data <- class.ind(vec)
					 			list_of_data[[i]] <- list_of_data[[i]][ ,!(names(list_of_data[[i]]) %in% covariates[i])]
					 			list_of_data[[i]] <- as.data.frame(cbind(list_of_data[[i]], new_data))
							}
					 	}
					 } else {
					 	for(j in 1:length(covariates)) {
					 			vec <- list_of_data[[i]][ ,covariates[j]]
					 			if(class(vec) %in% "character" | class(vec) %in% "factor") {
					 				new_vec <- as.numeric(vec)
					 				list_of_data[[i]] <- list_of_data[[i]][ ,!(names(list_of_data[[i]]) %in% covariates[i])]
					 				list_of_data[[i]] <- as.data.frame(cbind(list_of_data[[i]], new_vec))
					 			}
					 		}
					 }
					 object@list_of_data[[which_data[i]]] <- list_of_data[[i]]
				}
			} else {
				list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
				for(i in 1:length(which_data)) {
					dummy <- object@list_of_row_data[[which_data[i]]]
					dummy <- dummy[ ,covariates]
					dummy <- dummy[order(rownames(dummy)), ]
					data <- list_of_data[[i]]
					data <- data[rownames(data) %in% rownames(dummy), ]
					data <- data[order(rownames(data)), ]
					list_of_data[[i]] <- as.data.frame(cbind(data, dummy))
					if(class.ind) {
					 	for(j in 1:length(covariates)) {
					 		vec <- list_of_data[[i]][ ,covariates[j]]
					 		if(class(vec) %in% "character" | class(vec) %in% "factor") {
					 			new_data <- class.ind(vec)
					 			list_of_data[[i]] <- list_of_data[[i]][ ,!(names(list_of_data[[i]]) %in% covariates[i])]
					 			list_of_data[[i]] <- cbind.data.frame(list_of_data[[i]], new_data)
							}
					 	}
					 } else {
					 	for(j in 1:length(covariates)) {
					 			vec <- list_of_data[[i]][ ,covariates[j]]
					 			if(class(vec) %in% "character" | class(vec) %in% "factor") {
					 				new_vec <- as.numeric(vec)
					 				list_of_data[[i]] <- list_of_data[[i]][ ,!(names(list_of_data[[i]]) %in% covariates[i])]
					 				list_of_data[[i]] <- cbind.data.frame(list_of_data[[i]], new_vec)
					 			}
					 		}
					 }
					 object@list_of_data[[which_data[i]]] <- list_of_data[[i]]
				}
			}
			out <- object
			out <- object %>% add_function_info(function_name="add_phenotypes_as_covariates", 
				params=list(which_data=which_data, covariates=covariates, class.ind=class.ind, 
					phenotype=phenotype))
			return(out)
	})	




