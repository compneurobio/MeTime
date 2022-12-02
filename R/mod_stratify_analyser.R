
#' Function to stratify data in the metime analyser object
#' @description Function to stratify the data of interest into different objects that can be used
#' to perform calculations according the said stratification variable
#' @param object S4 object of class metime_analyser
#' @param which_data Dataset/datasets to be used for stratification
#' @param variable Phenotype based on which the stratification would be performed
#' @return list of metime_analyser objects which are stratified based on the variable chosen
#' @export
setGeneric("mod_stratify_analyser", function(object, which_data, variable) standardGeneric("mod_stratify_analyser"))
setMethod("mod_stratify_analyser", "metime_analyser", function(object, which_data, variable) {
				out <- lapply(which_data, function(x) {
								data <- object@list_of_data[[x]]
								row_data <- object@list_of_row_data[[x]]
								strat <- as.data.frame(nnet::class.ind(row_data[ ,variable]))
								rownames(strat) <- rownames(row_data)
								strat_list <- lapply(colnames(strat), function(y) {
												strat_samples <- ifelse(strat[,y]==1, rownames(strat), NA)
												strat_samples <- na.omit(strat_samples)
												strat_data <- data[rownames(data) %in% strat_samples, ]
												strat_row_data <- row_data[rownames(row_data) %in% strat_samples, ]
												object@list_of_data[[x]] <- strat_data
												object@list_of_row_data[[x]] <- strat_row_data
												other_data <- object@annotations
												other_object_data <- lapply(other_data, function(x) {
															new_data <- object@list_of_data[[x]]
															new_data <- new_data[rownames(new_data) %in% strat_samples, ]
															new_row_data <- object@list_of_row_data[[x]]
															new_row_data <- new_row_data[rownames(new_row_data) %in% strat_samples, ]
															object@list_of_data[[x]] <<- new_data
															object@list_of_row_data[[x]] <<- new_row_data
															return(NULL)
 													})
												return(object)
									})
								names(strat_list) <- paste(variable, colnames(strat), sep=": ")
								return(strat_list)
					})
				names(out) <- which_data
				return(out)
	})

