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


#' creating metime_plotter class that converts calculations and metadata as a plotable object to parse 
#' into viz_plotter
#' Contains slots - plot_data: Dataframe with plotting data and metadata for visualization
#' 				  - plot: ggplot(), circos() or visNetwork() object with predefined aesthetics 
#'                - calc_type: A vector to specify type of calculation - will be used for comp_ functions
#'                - calc_info: string to define the information about calculation
#' 				  - plot_type: A character vector to define the type of plots that are needed.
#' @rdname metime_plotter
#' @export
setClass("metime_plotter", slots=list(plot_data="list", plot="list", calc_type="character", calc_info="character", plot_type="character", style="character"))

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
		for(i in 1:length(object@list_of_data)) {
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
		}
	})


#' Setting new print definition for the metime_plotter object
#' @description function to see the structure of metime_plotter object
#' @param object S4 object of class metime_plotter
#' @examples structure(object)
#' @return structure of the S4 object
#' @export
setMethod("show", "metime_plotter", function(object) return(str(object, max.level=4)))

#' Setting new structure definition for the metime_plotter object
#' @description function to see the structure of metime_plotter object
#' @param object S4 object of class metime_plotter
#' @examples structure(object)
#' @return structure of the S4 object
#' @export
setMethod("structure", "metime_plotter", function(object) return(str(object, max.level = 3)))

