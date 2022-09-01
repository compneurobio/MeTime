usethis::use_package("tidyverse", type="depends")
usethis::use_package("igraph", type="depends")
usethis::use_package("plyr", type="depends")
usethis::use_package("DT", type="depends")
usethis::use_package("M3C", type="depends")
usethis::use_package("umap", type="depends")
usethis::use_package("ggplot2", type="depends")
usethis::use_package("plotly", type="depends")
usethis::use_package("RColorBrewer", type="depends")
usethis::use_package("GeneNet", type="depends")
usethis::use_package("longitudinal", type="depends")
usethis::use_package("glmnet", type="depends")
usethis::use_package("Boruta", type="depends")
usethis::use_package("nnet", type="depends")
usethis::use_package("shiny", type="depends")
usethis::use_package("abind", type="depends")
usethis::use_package("multiway", type="depends")

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
setMethod("show", "metime_analyser", function(object) return(str(object, max.level=3)))


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
