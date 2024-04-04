#' Get metadata for columns(in most cases for metabolites)
#' @description function to generate a metadata list for building the MeTime plotter object
#' @param object a S4 object of the class "metime_analyzer".
#' @param which_data a character to define which dataset is to be used.
#' @param columns A named list of a vector of named characters containing the columns of interest. Length of the list should be same as length of which_data. List can be named or unnamed but the character vectors should be named.
#' Moreover, they should have same names. Default is set to NULL which results in empty metadata dataframe.
#' Ex: list(nmr_data=c(id="id", sub_pathway="Group"), lipid_data=c(id="id", sub_pathway="sub_pathway"))
#' @return a dataframe with metadata information.
#' @export
setGeneric("get_metadata_for_columns", function(object, which_data, columns=NULL) standardGeneric("get_metadata_for_columns"))
setMethod("get_metadata_for_columns", "metime_analyser", function(object, which_data, columns=NULL) {
				stopifnot(which_data %in% names(object@list_of_data))
				
				if(is.null(columns)) return(NULL)

				names <- lapply(seq_along(columns), function(x) return(names(columns[[x]])))
				names_check <- Reduce(intersect, names)
				
				lapply(seq_along(names), function(i) {
					if(!all(names[[i]] %in% names_check)) {
						stop("names of the character vectors in columns do not match")
					}
				})

				col_data <- lapply(which_data, function(a) {
							return(get_coldata(object=object, which_data=a))
						})
				

				if(!"id" %in% names_check) {
					out <- lapply(seq_along(col_data), function(a) {
							data <- col_data[[a]][ ,c("id", columns[[a]])]
							colnames(data) <- c("id", names_check)
							rownames(data) <- data$id
							data$class <- rep(which_data[a], each=length(rownames(data)))
							data[order(rownames(data)),]
							return(data)
						}) %>% do.call(what=rbind.data.frame)
				} else {
					out <- lapply(seq_along(col_data), function(a) {
							data <- col_data[[a]][ ,columns[[a]]]
							colnames(data) <- names_check
							rownames(data) <- data$id
							data$class <- rep(which_data[a], each=length(rownames(data)))
							data <- data[order(rownames(data)), ]
							return(data)
						}) %>% do.call(what=rbind.data.frame)
				}
				return(out)
	})

#metadata <- get_metadata_for_columns(object=object, which_data="lipid_data", columns=list(c("id", "sub_pathway")), 
#									names=c("name", "group"), index_of_names="id")
