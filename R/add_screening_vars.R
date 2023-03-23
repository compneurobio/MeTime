#' Add values measured at a screening time for samples to be added to all time points
#' @description Add all data points that were measured only during a defined screening time point to all other measurements of the same subject.
#' Note that this function is only valid in the case of phenotype data and not in the case of row_data
#' @param object a S4 object of class "metime_analyser".
#' @param vars a character vector naming the columns of interest
#' @return S4 object of class "metime_analyser" with screening measurements added to all other time points of the same subject.
#' @export
setGeneric("add_screening_vars", function(object, vars) standardGeneric("add_screening_vars"))

setMethod("add_screening_vars", "metime_analyser", function(object, vars) {
	phenotype_name <- object@annotations[[1]]$phenotype
	phenotype <- object@list_of_data[[phenotype_name]]
	if(!all(vars %in% names(phenotype))) {
		warning("All variables mentioned are not in the phenotype dataset. Exiting without making any changes")
		return(object)
	}
	phenotype$subject <- rownames(phenotype) %>% gsub(pattern="_[a-z|A-Z][-|0-9]+", replacement="")
	screening <- phenotype[grep("-1|-2", rownames(phenotype)), ]
	new_rows <- as.data.frame(screening[, c(vars, "subject")])
	new_rows <- na.omit(new_rows)
	new_rows <- new_rows[order(rownames(new_rows)), ]
	sample_names <- new_rows$subject
	new_phenotype <- lapply(seq_along(sample_names), function(i) {
		new_data <- phenotype[phenotype$subject %in% sample_names[i], ] %>% as.data.frame()
		new_data[ ,vars] <- lapply(new_rows[new_rows$subject %in% sample_names[i], vars], 
						rep, length(new_data$subject)) %>% 
						as.data.frame()
		return(new_data)
		}) %>% do.call(what=rbind.data.frame)
	new_phenotype <- new_phenotype[order(rownames(new_phenotype)), ] 
	object@list_of_data[[phenotype_name]] <- new_phenotype
	out <- object
	out <- add_function_info(object=out, function_name="add_screening_vars", params=list(vars=vars))
	return(out)
})




