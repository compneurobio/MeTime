#' Function to add measurements taken at screening time for samples to be added to all timepoints
#' @description A method applied on the s4 object of class "metime_analyser" to add all those datapoints that were measured only during screening
#' to all the respective samples at all timepoints
#' @examples # adding APOEGrp, PTGENDER to all data points
#' new_with_apoegrp_sex <- add_screening_vars(object=metime_analyser_object, vars=c("APOEGrp","PTGENDER"))
#' @param object An object of class metime_analyser
#' @param vars A character naming the vars of interest
#' @return phenotype data which can be replaced into the original object or use it separately with a different object
#' @export
setGeneric("add_screening_vars", function(object, vars) standardGeneric("add_screening_vars"))

setMethod("add_screening_vars", "metime_analyser", function(object, vars) {
	phenotype_name <- object@annotations[[1]]$phenotype
	phenotype <- object@list_of_data[[phenotype_name]]
	screening <- phenotype[grep("-1|-2", rownames(phenotype)), ]
	new_rows <- as.data.frame(screening[, vars])
	new_rows <- na.omit(new_rows)
	new_rows <- new_rows[order(rownames(new_rows)), ]
	sample_names <- rownames(new_rows) %>% gsub(pattern="_[a-z|A-Z][0-9]+", replacement="")
	test <- lapply(seq_along(sample_names), function(i) {
		samples <- rownames(phenotype) %>% gsub(pattern="_[a-z|A-Z][0-9]+", replacement="")
		index <- which(samples %in% sample_names[i], arr.ind=TRUE)
		phenotype <- lapply(index, function(j) {
				phenotype[j, vars] <- new_rows[rownames(new_rows)[i], vars]
			}) 
		})
	object@list_of_data[[phenotype_name]] <- phenotype
	out <- object
	out <- add_function_info(object=out, function_name="add_screening_vars", params=list(vars=vars))
	return(out)
})




