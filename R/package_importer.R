usethis::use_import_from("rstatix", fun="t_test")
usethis::use_package("tidyverse", type="depends")
usethis::use_package("igraph", type="Imports")
usethis::use_package("plyr", type="Imports")
usethis::use_package("DT", type="Imports")
usethis::use_import_from("M3C", fun="tsne")
usethis::use_import_from("umap", fun="umap")
usethis::use_package("plotly", type="Imports")
usethis::use_package("visNetwork", type="Imports")
usethis::use_package("RColorBrewer", type="Imports")
usethis::use_package("GeneNet", type="Imports")
usethis::use_package("longitudinal", type="Imports")
usethis::use_package("glmnet", type="Imports")
usethis::use_package("Boruta", type="Imports")
usethis::use_import_from("nnet", fun="class.ind")
usethis::use_import_from("abind", fun="abind")
usethis::use_import_from("multiway", fun="parafac")
usethis::use_import_from("DT", fun="datatable")
usethis::use_import_from("xlsx", fun="write.xlsx")
usethis::use_package("visNetwork", type="Imports")
usethis::use_package("mgcv", type="Imports")
usethis::use_import_from("lmerTest", fun="lmer")
usethis::use_import_from("parallel", fun="mclapply")
usethis::use_import_from("reshape2", fun="melt")
usethis::use_import_from("data.table", fun="setDT")
usethis::use_package("MatrixEQTL", type="Imports")
usethis::use_package("WGCNA", type="Imports")
usethis::use_import_from("dynamicTreeCut", fun="cutreeDynamic")


#' Validity function to check if the object is valid or not
#' @description Function to check the validity of metime_analyser 
#' @param object An S4 object of class metime_analyser
#' @examples validObject(object)
#' @returns logical suggesting if the object is intact or not
#' @export
validity <- function(object) {
		out <- TRUE
		if(!check_rownames_and_columns(object)) out <- FALSE 
		if(!check_ids_and_classes(object)) out <- FALSE
		if(!check_results(object)) out <- FALSE
		return(out) 
	}

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
								 annotations="list", results="list"), validity=validity) 



# add a function to generate Rmarkdown directly from the analyser object - viz_analyser()
# create another function which add the function applied on the analyser_object - ongoing
# Try to improve the code of all functions and then Matthias will check it

#' Setting a plotting method for the metime_analyser class
#' @description Function to plot results of a certain calculation 
#' @param object An S4 object of class metime_analyser
#' @param results_index Index/name of the results to be plotted
#' @return plots for a certain set of results
#' @export
setGeneric("plot", function(object, results_index, ...) standardGeneric("plot"))
setMethod("plot", "metime_analyser", function(object, results_index, ...) {
			results <- object@results[[results_index]]
	})


#' Setting new structure definition for the metime_analyser object
#' @description function to see the structure of metime_analyser object
#' @param object S4 object of class metime_analyser
#' @examples structure(object)
#' @return structure of the S4 object
#' @export
setGeneric("structure", function(object) standardGeneric("structure"))
setMethod("structure", "metime_analyser", function(object) return(str(object, max.level = 3)))

#' Setting new print definition for the metime_analyser object
#' @description function to see the structure of metime_analyser object
#' @param object S4 object of class metime_analyser
#' @examples structure(object)
#' @return structure of the S4 object
#' @export
setMethod("show", "metime_analyser", function(object) {
		out <- lapply(seq_along(object@list_of_data), function(i) {
			if(names(object@list_of_data)[i] %in% object@annotations$phenotype) {
				cat(paste("Dataset:", names(object@list_of_data)[i], "with", 
				dim(object@list_of_data[[i]])[1], "samples", "and", dim(object@list_of_data[[i]])[2], "phenotypes",sep=" "))
				cat("\n")
			} else if(names(object@list_of_data)[i] %in% object@annotations$medication) {
				cat(paste("Dataset:", names(object@list_of_data)[i], "with", 
				dim(object@list_of_data[[i]])[1], "samples", "and", dim(object@list_of_data[[i]])[2], "medicines", sep=" "))
				cat("\n")
			} else {
				cat(paste("Dataset:", names(object@list_of_data)[i], "with", 
				dim(object@list_of_data[[i]])[1], "samples", "and", dim(object@list_of_data[[i]])[2], "metabolites",sep=" "))
				cat("\n")
			}
			return(NULL)
		})
	})
