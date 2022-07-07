
#creating reference metab-analyser class that creates an object with full data

#' Constructor to generate an object of class metab_analyser. 
#' contains slots - list_of_data: For the list of all data matrices.
#'				  - list_of_col_data: list of all the col data files in the same order.
#' 				  - list_of_row_data: list of all the row data files in the same order.
#' 				  - annotations: list with phenotype and medication. Each of which is character that represents 
#'									the name of the aforementioned dataset types.  	
#' 	
#' @rdname metab_analyser
#' @export 
setClass("metab_analyser", slots=list(list_of_data="list", list_of_col_data="list", list_of_row_data="list", 
								 annotations="list")) 

#' Function to add measurements taken at screening time for samples to be added to all timepoints
#' @description A method applied on the s4 object of class "metab_analyser" to add all those datapoints that were measured only during screening
#' to all the respective samples at all timepoints
#' @param object An object of class metab_analyser
#' @param vars A character naming the vars of interest
#' @return phenotype data which can be replaced into the original object or use it separately with a different object
#' @export
setGeneric("add_screening_vars", function(object, vars) standardGeneric("add_screening_vars"))

setMethod("add_screening_vars", "metab_analyser", function(object, vars) {
		phenotype_name <- object@annotations$phenotype
		phenotype <- object@list_of_data[[phenotype_name]]
		screening <- phenotype[grep("-1|-2", rownames(phenotype)), ]
		new_rows <- as.data.frame(screening[, vars])
		new_rows <- na.omit(new_rows)
		sample_names <- unlist(lapply(strsplit(rownames(screening), split="_t"), function(x) return(x[1])))
		for(i in 1:length(sample_names)) {
			index <- grep(sample_names[i], rownames(phenotype)) 
			for(j in index) {
				phenotype[j, vars] <- screening[rownames(screening)[i], vars]
			}
		}
		return(phenotype)
	})

