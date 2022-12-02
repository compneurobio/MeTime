#' Function to add measurements taken at screening time for samples to be added to all timepoints in row data
#' @description A method applied on the s4 object of class "metime_analyser" to add all those datapoints that were measured only during screening
#' to all the respective samples at all timepoints in row_data lists
#' @examples # adding APOEGrp, PTGENDER, and diag group to all data points and prepping the object for viz_distribution_plotter()
#' object <- add_distribution_vars_to_rows(object=data, screening_vars=c("APOEGrp", "DXGrp_longi", "PTGENDER"), 
#'			distribution_vars=c("Age", "BMI", "ADNI_MEM", "ADNI_LAN", "ADNI_EF", "APOEGrp", "DXGrp_longi", "PTGENDER"), which_data="lipid_data")
#' @param object An object of class metime_analyser
#' @param vars A character naming the vars of interest
#' @param which_data dataset to which the information is to be added(only 1 can be used at a time)
#' @return object of class metime_analyser with phenotype data added to row data
#' @export
setGeneric("add_distribution_vars_to_rows", function(object, screening_vars, distribution_vars, which_data) standardGeneric("add_distribution_vars_to_rows"))
setMethod("add_distribution_vars_to_rows", "metime_analyser", function(object, screening_vars, distribution_vars, which_data) {
			if(!is.null(screening_vars)) {
				phenotype <- add_screening_vars(object, screening_vars)
			} else {
				phenotype <- object@list_of_data[[object@annotations$phenotype]]
			}
			data <- as.data.frame(object@list_of_row_data[[which_data]])
			phenotype <- phenotype[rownames(phenotype) %in% rownames(data), ]
			phen_data <- phenotype[order(rownames(phenotype)),]
			data <- data[order(rownames(data)), ]
			data <- cbind(data, phen_data[ ,distribution_vars])
			object@list_of_row_data[[which_data]] <- data
			out <- object
			return(out)
	})

