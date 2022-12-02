
#' Function to convert metabolite names to IDs
#' @description Function to convert metabolite names to IDs 
#' @param object An S4 object of class metime_analyser
#' @param which_data character vector to define the datasets to use
#' @return A list with S4 object and list of mapping tables, the object can be used for GGMs
#' @export
setGeneric("mod_code_metab_names", function(object, which_data) standardGeneric("mod_code_metab_names"))
setMethod("mod_code_metab_names", "metime_analyser", function(object, which_data) {
					tables <- list()	
					for(i in 1:length(which_data)) {
								data <- object@list_of_data[[which_data[i]]]
								metabs <- paste(unlist(strsplit(which_data[i], split=""))[1], 1:length(colnames(data)), sep=".")
								tables[[i]] <- as.data.frame(cbind(colnames(data), metabs))
								colnames(tables[[i]]) <- c("id", "metabolite")
								colnames(data) <- metabs
								object@list_of_data[[which_data[i]]] <- data
					}
					table <- as.data.frame(do.call(rbind, tables))
					out <- list(object=object, table=table)
					return(out)
		})


