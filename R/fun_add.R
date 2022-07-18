
#creating reference metime-analyser class that creates an object with full data

#' Constructor to generate an object of class metime_analyser. 
#' contains slots - list_of_data: For the list of all data matrices.
#'				  - list_of_col_data: list of all the col data files in the same order.
#' 				  - list_of_row_data: list of all the row data files in the same order.
#' 				  - annotations: list with phenotype and medication. Each of which is character that represents 
#'									the name of the aforementioned dataset types.  	
#' 	
#' @rdname metime_analyser
#' @export 
setClass("metime_analyser", slots=list(list_of_data="list", list_of_col_data="list", list_of_row_data="list", 
								 annotations="list")) 

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

#' Function to check normality and add data to col data
#' @description A method applied on the s4 object of class "metime_analyser" to check normality of the metabolites
#' and add it to corresponding columns
#' @examples 
#'	object <- add_col_normality(object=data, which_data=c("lipid_data","nmr_data"), type="shapiro", metab_names=c("metabolite","id"))
#' @param object An object of class metime_analyser
#' @param which_data dataset on which the method is to be applied
#' @param type type of test, currently only shapiro-wilk test is available under the choice "shapiro"
#' @param metab_names column that has the metabolite names in col_data.
#' @return S4 object with shapiro wilk test related data in the col_data
#' @export
setGeneric("add_col_normality", function(object, which_data, type, metab_names) standardGeneric("add_col_normality"))
setMethod("add_col_normality", "metime_analyser", function(object, which_data, type="shapiro", metab_names) {
	if(type %in% "shapiro") {
      	for(i in 1:length(which_data)) {
      		out <- lapply(names(object@list_of_data[[which_data[i]]]), function(x) {
      			my_model <- shapiro.test(object@list_of_data[[which_data[i]]][[x]])
      			return(data.frame(id=x,
                 	shapiro_pval=as.numeric(my_model$p.value),
                 	shapiro_statistic=as.numeric(my_model$statistic),
                 	stringsAsFactors = F))}) %>%  
      				do.call(what=rbind.data.frame)
      		dummy <- out[order(out$id), ]
      		shapiro_pval <- dummy$shapiro_pval
      		shapiro_normal <- dummy$shapiro_normal
      		shapiro_statistic <- dummy$shapiro_statistic
      		shapiro_normal <- ifelse(dummy$shapiro_pval>0.05, TRUE,FALSE)
      		col_data <- object@list_of_col_data[[which_data[i]]]
      		col_data <- col_data[order(col_data[ ,metab_names[i]]), ]
      		col_data <- as.data.frame(cbind(col_data, shapiro_pval, shapiro_statistic, shapiro_normal))
      		object@list_of_col_data[[which_data[i]]] <- col_data  
      	}
      	
	}
	return(object)
  
})

