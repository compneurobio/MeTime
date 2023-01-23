
#' Function to check the ids in the data and data format
#' @description sanity check to check for ids and order of the data
#' @param object S4 object of class of metime_analyser
#' @return NULL if it passes all the sanity checks
#' @export
setGeneric("check_ids_and_classes",  function(object) standardGeneric("check_ids_and_classes")) 
setMethod("check_ids_and_classes", "metime_analyser", function(object) {
			list_of_interest <- list(data=object@list_of_data, 
									col_data=object@list_of_col_data, 
									row_data=object@list_of_row_data)
			names <- names(object@list_of_data)
			type <- unlist(lapply(strsplit(names(list_of_interest), split="list_of_"), function(x) return(x[2])))
			out <- TRUE
			test <- lapply(seq_along(list_of_interest), function(x) {
						lapply(seq_along(list_of_interest[[x]]), function(y) {
							tryCatch(class(list_of_interest[[x]][[y]]) %in% "data.frame", 
								error=function(e) {
									message("In dataset ", names[y], " datatype ", type[x], 
										" class of data should be data.frame only")
									out <- FALSE
								})
							tryCatch(all(rownames(list_of_interest[[x]][[y]]) == list_of_interest[[x]][[y]]$id), 
								error=function(e) {
									message("In dataset", names[y], "datatype", type[x], 
										"rownames and ids are not mapped correctly")
									out <- FALSE
								})
						})
					})
			return(out)
	})

#' Function to check the format of rownames and colnames and if they are same or not
#' @description sanity check to check for rownames of the data
#' @param object S4 object of class of metime_analyser
#' @return NULL if it passes all the sanity checks
#' @export
setGeneric("check_rownames_and_colnames", function(object) standardGeneric("check_rownames_and_colnames"))
setMethod("check_rownames_and_colnames", "metime_analyser", function(object) {
		list_of_data <- object@list_of_data
		list_of_row_data <- object@list_of_row_data
		list_of_col_data <- object@list_of_col_data
		out <- TRUE
		#data_checks for duplicates
		test <- lapply(seq_along(list_of_data), function(i) {
					tryCatch(length(colnames(list_of_data[[i]])) == length(unique(colnames(list_of_data[[i]]))),
						error=function(e) {
							message("In Dataset ", names(list_of_data)[i], " some metabolites are duplicated", sep=" ")
							out <- FALSE
						})
					tryCatch(length(rownames(list_of_data[[i]])) == length(unique(rownames(list_of_data[[i]]))), 
						error=function(e) {
							message("In Dataset ", names(list_of_data)[i], " some samples are duplicated", sep=" ")
							out <- FALSE
						})
					if(!(names(list_of_data)[i] %in% object@annotations$phenotype | names(list_of_data)[i] %in% object@annotations$medication)) {
						tryCatch(all(list_of_col_data[[i]]$id %in% colnames(list_of_data[[i]])), 
							error=function(e) {
								message("In Dataset ", names(list_of_data)[i], " all metabolites are not mapped in the col data")
								out <- FALSE
							})
						tryCatch(all(list_of_row_data[[i]]$id %in% rownames(list_of_data[[i]])), 
							error=function(e) {
								out <- FALSE
								message("In Dataset ", names(list_of_data)[i], " all samples are not mapped in the row data")
							})
					}
					return(NULL)
				})
		return(out)
	})


#' Function to check the format of results if they exist
#' @description sanity check to check for results of the analysis
#' @param object S4 object of class of metime_analyser
#' @return NULL if it passes all the sanity checks
#' @export
setGeneric("check_results", function(object) standardGeneric("check_results"))
setMethod("check_results", "metime_analyser", function(object) {
		out <- TRUE
		if(!length(object@results)==0) {
			test <- lapply(seq_along(object@results), function(x) {
					if(length(grep("calc_", object@results[[x]][["functions"]]))>=1) {
						tryCatch(object@results[[x]][["plot"]] %>% is_ggplot() |
							all(class(object@results[[x]][["plot"]]) %in% c("plotly", "htmlwidget")) |
							  all(class(object@results[[x]][["plot"]]) %in% c("visNetwork", "htmlwidget")), 
							  	error=function(e) {
							  			message("plot of the results in", names(object@results)[x], " is missing")
							  			out <<- FALSE
							  		})
						if(class(object@results[[x]][["plot_data"]]) %in% "data.frame") {
							tryCatch(dim(object@results[[x]][["plot_data"]])[1] != 0, 
								error= function(e) {
									message("plot_data is wrong in ", names(object@results)[x])
								})
						} else if(class(object@results[[x]][["plot_data"]]) %in% "list") {
							tryCatch(dim(object@results[[x]][["plot_data"]][["node"]])[1] != 0 && 
										dim(object@results[[x]][["plot_data"]][["edge"]])[1] != 0, 
								error= function(e) {
									message("plot_data is wrong in ", names(object@results)[x])
								})
						}
					}
				})
		}
		return(out) 
	})