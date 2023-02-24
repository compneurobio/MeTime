
#' Function for Plotting distributions of phenotypic variables 
#' @description A method to be applied onto s4 object so as to obtain distributions of various phenotypic variables
#' @param object An object of class metime_analyser
#' @param cols character vector to define the columns whose distributions are wanted from row_data
#' @param which_data Name of the dataset from which the samples will be extracted
#' @param stratifications List to define the stratification of interest
#' @param name Character to name the results 
#' @return S4 object with updated plot_data and plots 
#' with plots being either 1) density plot
#'							2) bar plot and a line plot
#' @export
setGeneric("calc_distribution_samples", function(object, which_data, cols, stratifications, name="calc_distribution_samples_1") standardGeneric("calc_distribution_samples"))
setMethod("calc_distribution_samples", "metime_analyser",function(object, which_data, cols, stratifications=NULL, name="calc_distribution_samples_1") {
	
	if(grep(name, names(object@results)) %>% length() >=1) {
    	warning("name of the results was previously used, using a different name")
    	index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
    	index <- c(0:9)[grep(index, 0:9)+1]
    	name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
  	}
	stopifnot(length(which_data)==1)
	
	if(!is.null(stratifications)) {
		data_list <- get_stratified_data(object=object, 
			which_data=which_data, stratifications=stratifications)
    	row_data <- data_list[["row_data"]]
	} else {
		row_data <- get_rowdata(object, which_data=which_data)
	}
	
    
    if(is.null(cols)) {
    	row_data <- row_data[ ,c("subject", "time")]
    } else {
    	row_data <- row_data[ ,c(cols, "subject", "time")]
    }
    #### results data frames to be built here
    out <- get_make_results(object=object, data=list(row_data), 
    	metadata=NULL, calc_type="distribution_samples", 
    	calc_info=paste("Samples distributions dataframe of columns of interest in", which_data, sep=" "), 
    	name=name) %>%
    	add_function_info(function_name="calc_distribution_samples", 
    		params=list(which_data=which_data, cols_for_meta=cols,
    			stratifications=stratifications, name=name)) %>% update_plots(type="distribution")
    return(out)
})



#### SPLIT Calc_distibution and calc_dimensionality_reduction
#' Function for Plotting distributions of phenotypic variables 
#' @description A method to be applied onto s4 object so as to obtain distributions of various phenotypic variables
#' @param object An object of class metime_analyser
#' @param which_data A Character vector to define the dataset of interest
#' @param cols character vector to define the columns whose distributions are wanted from col_data#' @param which_data Name of the dataset from which the samples will be extracted
#' @param name Character to name the results 
#' @return S4 object with updated plot_data and plots 
#' with plots being either 1) density plot
#'							2) bar plot and a line plot
#' @export
setGeneric("calc_distribution_metabs", function(object, which_data, cols, name="calc_distribution_metabs_1") standardGeneric("calc_distribution_metabs"))
setMethod("calc_distribution_metabs", "metime_analyser",function(object, which_data, cols, name="calc_distribution_metabs_1") {
	if(grep(name, names(object@results)) %>% length() >=1) {
    	warning("name of the results was previously used, using a different name")
    	index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
    	index <- c(0:9)[grep(index, 0:9)+1]
    	name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
  	}
	stopifnot(length(which_data)==1)
	col_data <- get_coldata(object=object, which_data=which_data)
	if(is.null(cols)) {
    	col_data <- col_data$id
    } else {
    	col_data <- col_data[ ,c(cols, "id")]
    }
    out <- get_make_results(object=object, data=list(col_data), 
    	metadata=NULL, calc_type="distribution_metabs", 
    	calc_info=paste("Metabolites distributions dataframe of columns of interest in", which_data, sep=" "), 
    	name=name) %>%
    	add_function_info(function_name="calc_distribution_metabs", 
    		params=list(which_data=which_data, cols_for_meta=cols,
    			name=name)) %>% update_plots(type="distribution")
})

