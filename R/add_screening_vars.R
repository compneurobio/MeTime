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
	phenotype$subject <- rownames(phenotype) %>% gsub(pattern="_[a-z|A-Z][-|0-9]+", replacement="")
	screening <- phenotype[grep("-1|-2", rownames(phenotype)), ]
	new_rows <- as.data.frame(screening[, c(vars, "subject")])
	new_rows <- na.omit(new_rows)
	new_rows <- new_rows[order(rownames(new_rows)), ]
	sample_names <- new_rows$subject
	new_phenotype <- lapply(seq_along(sample_names), function(i) {
		samples <- phenotype$subject
		new_data <- phenotype[sample_names[i] %in% phenotype$subject, ] %>% as.data.frame()
		new_data[ ,vars] <- lapply(new_rows[new_rows$subject %in% sample_names[i], vars], rep, length(new_data$subject)) %>% 
						as.data.frame() %>% return()
		}) %>% do.call(what=rbind.data.frame)
	new_phenotype <- new_phenotype[order(rownames(new_phenotype)), ] 
	object@list_of_data[[phenotype_name]] <- new_phenotype
	out <- object
	out <- add_function_info(object=out, function_name="add_screening_vars", params=list(vars=vars))
	return(out)
})




