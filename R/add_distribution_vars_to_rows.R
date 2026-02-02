#' Add phenotypic measurements that are not added to the row_data of the dataset
#' @description A method applied on the s4 object of class "metime_analyser" to add all those datapoints that are present in phenotype dataframe but not in row_data
#' comes with the feature of updating those data points measured only at screening to all datapoints and then adding it to row_data
#' @examples 
#' object <- add_distribution_vars_to_rows(object=data, screening_vars=c("var1", "var2"), 
#'			distribution_vars=c("var1", "var2", "var3"), which_data="dataset1")
#' @param object An object of class metime_analyser
#' @param distribution_vars A character naming the vars of interest
#' @param screening_vars A character vector to define the vars that are to be updated as per add_screening_vars() else set it to NULL(Default).
#' @param which_data dataset to which the information is to be added(only 1 can be used at a time)
#' @seealso [add_screening_vars]
#' @return object of class metime_analyser with phenotype data added to row data
#' @export
setGeneric("add_distribution_vars_to_rows", function(object, screening_vars, distribution_vars, which_data) standardGeneric("add_distribution_vars_to_rows"))
setMethod("add_distribution_vars_to_rows", "metime_analyser", function(object, screening_vars, distribution_vars, which_data) {
			if(!which_data %in% names(object@list_of_data)) {
				warning("dataset not found in the object. Exiting without making any changes")
				return(object)
			}
			if(!is.null(screening_vars)) {
				object <- add_screening_vars(object, screening_vars)
			} 
			phenotype <- object@list_of_data[[object@annotations[[1]]$phenotype]]
			if(!all(distribution_vars %in% names(phenotype))) {
				warning("All variables mentioned are not in the phenotype dataset. Exiting without making any changes")
				return(object)
			}
			data <- as.data.frame(object@list_of_row_data[[which_data]])
			phenotype <- phenotype[rownames(phenotype) %in% rownames(data), ]
			phen_data <- phenotype[order(rownames(phenotype)),]
			data <- data[order(rownames(data)), ]
			data <- cbind.data.frame(data, phen_data[ ,distribution_vars])
			object@list_of_row_data[[which_data]] <- data
			out <- object
			out <- add_function_info(object=out, function_name="add_distribution_vars_to_rows",
					params=list(screening_vars=screening_vars, distribution_vars=distribution_vars, which_data=which_data))
			return(out)
	})

