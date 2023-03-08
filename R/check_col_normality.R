
#' Function to check for col_normality data whether it is added or not. 
#' @description function to check whether col_normality data is added to
#' the object or not
#' @param object S4 object of class of metime_analyser
#' @param which_data dataset/s to check
#' @importClassesFrom metime_analyser
#' @return NULL if it passes all the sanity checks
#' @export
setGeneric("check_col_normality", function(object, which_data) standardGeneric("check_col_normality"))
setMethod("check_col_normality", "metime_analyser", function(object, which_data) {
			check_data <- object@list_of_col_data[names(object@list_of_col_data) %in% which_data]
			cols <- c("pval", "statistic", "normal")
			test <- lapply(seq_along(check_data), function(i) {
						if(is.null(grep(paste(cols, collapse="|"), colnames(check_data[[i]])))) {
							print(paste("normality data is not added to", which_data[i], "dataset", sep=" "))
						} else {
							if(is.null(grep("shapiro", colnames(check_data[[i]])))) {
								print(paste("shapiro normality test is not added to", which_data[i], "dataset", sep=" "))
							} else if(is.null(grep("kruskal", colnames(check_data[[i]])))) {
								print(paste("kruskal normality test is not added to", which_data[i], "dataset", sep=" "))
							} else {
								print(paste("All normality tests are added to", which_data[i], "dataset", sep=" "))
							}
						}
				})
	})

