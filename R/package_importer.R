usethis::use_import_from("rstatix", fun="t_test")
usethis::use_package("tidyverse", type="depends")
usethis::use_package("igraph", type="Imports")
usethis::use_package("plyr", type="Imports")
usethis::use_package("DT", type="Imports")
usethis::use_import_from("M3C", fun="tsne")
usethis::use_import_from("umap", fun="umap")
usethis::use_package("ggplot2", type="Imports")
usethis::use_package("plotly", type="Imports")
usethis::use_package("visNetwork", type="Imports")
usethis::use_package("RColorBrewer", type="Imports")
usethis::use_package("GeneNet", type="Imports")
usethis::use_package("longitudinal", type="Imports")
usethis::use_package("glmnet", type="Imports")
usethis::use_package("Boruta", type="Imports")
usethis::use_import_from("nnet", fun="class.ind")
usethis::use_package("shiny", type="Imports")
usethis::use_import_from("abind", fun="abind")
usethis::use_import_from("multiway", fun="parafac")
usethis::use_import_from("DT", fun="datatable")
usethis::use_import_from("xlsx", fun="write.xlsx")
usethis::use_package("visNetwork", type="Imports")
usethis::use_package("mgcv", type="Imports")
usethis::use_package("lme4", type="Imports")

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

