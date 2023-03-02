


#' Get metadata for columns(in most cases for metabolites)
#' @description function to generate a metadata list for building the MeTime plotter object
#' @param object S4 object of class MeTime Analyser
#' @param which_data Names of dataset/s to be used
#' @param columns A list of character vectors for the columns of interest. Length of the list should be
#' same as length of which_data. List can be named or unnamed but the character vectors should be named.
#' Moreover, they should have same names. Default is set to NULL which results in empty metadata dataframe.
#' Ex: list(nmr_data=c(id="id", sub_pathway="Group"), lipid_data=c(id="id", sub_pathway="sub_pathway"))
#' @return data.frame with metadata information
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
							data$id <- NULL 
							data$class <- rep(which_data[a], each=length(data$id))
							return(data)
						}) %>% do.call(what=rbind.data.frame)
				} else {
					out <- lapply(seq_along(col_data), function(a) {
							data <- col_data[[a]][ ,columns[[a]]]
							colnames(data) <- names_check
							col_of_int <- names(columns[[a]])[which(columns[[a]] %in% "id")]
							rownames(data) <- data[ ,col_of_int]
							data[ ,col_of_int] <- NULL 
							data$class <- rep(which_data[a], each=length(data$id))
							return(data)
						}) %>% do.call(what=rbind.data.frame)
				}
				return(out)
	})

#metadata <- get_metadata_for_columns(object=object, which_data="lipid_data", columns=list(c("id", "sub_pathway")), 
#									names=c("name", "group"), index_of_names="id")


# list_of_col_data <- object@list_of_col_data[names(object@list_of_col_data) %in% which_data]
# 				if(length(which_data)>1) {
# 					if(is.null(columns)) {
# 						out <- NULL
# 						return(out)
# 					} else {
# 						list_of_metadata_metabs <- list()
# 						for(i in 1:length(columns)) {
# 							list_of_metadata_metabs[[i]] <- as.matrix(list_of_col_data[[i]][, columns[[i]]])
# 							class <- rep(which_data[i], each=length(list_of_metadata_metabs[[i]][,1]))
# 							list_of_metadata_metabs[[i]] <- list_of_metadata_metabs[[i]][order(list_of_metadata_metabs[[i]][ ,index_of_names[i]]), ]
# 							colnames(list_of_metadata_metabs[[i]]) <- names
# 							list_of_metadata_metabs[[i]] <- as.data.frame(cbind(list_of_metadata_metabs[[i]], class))
# 						}
# 						out <- lapply(list_of_metadata_metabs, as.data.frame)
# 						out <- as.data.frame(do.call(rbind, out))
# 						rownames(out) <- out[, names[1]] 
# 					}
# 				} else {
# 					if(is.null(columns)) {
# 						out <- NULL
# 						return(out)
# 					} else {
# 						col_data <- list_of_col_data[[1]]
# 						class <- rep(which_data, each=length(col_data[,1]))
# 						out <- col_data[ ,columns[[1]]]
# 						out <- as.data.frame(cbind(out, class))
# 						out <- out[order(out[,index_of_names]), ]
# 						rownames(out) <- out[,index_of_names]
# 						colnames(out) <- c(names, "class")
# 					}
# 				}