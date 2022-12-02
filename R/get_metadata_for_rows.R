#' Get metadata for rows(in most cases for samples)
#' @description function to generate a metadata list for building the MeTime plotter object
#' @param object S4 object of class MeTime Analyser
#' @param which_data Names of dataset/s to be used
#' @param columns A list of character vectors for the columns of interest. Length of the list should be
#' same as length of which_data
#' @return data.frame with metadata information for rows
#' @export
setGeneric("get_metadata_for_rows", function(object, which_data, columns) standardGeneric("get_metadata_for_rows"))
setMethod("get_metadata_for_rows", "metime_analyser", function(object, which_data, columns) {
					if(length(which_data) > 1) {
							if(is.null(columns)) {
									out <- NULL
									return(out)
							} else {
									object <- mod_extract_common_samples(object)
									list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
									list_of_data <- lapply(list_of_data, function(x) return(x[order(rownames(x)), ]))
									out <- object@list_of_row_data[[1]][ ,columns]
									timepoints <- unlist(lapply(strsplit(rownames(out), split="_"), function(x) return(x[2])))
									samples <- unlist(lapply(strsplit(rownames(out), split="_"), function(x) return(x[1])))
									levels <- sort(unique(as.numeric(unlist(lapply(strsplit(timepoints, split="t"), function(x) return(x[2]))))))
									timepoints <- factor(timepoints, levels=paste("t",levels,sep=""))
									out <- as.data.frame(cbind(out, timepoints, samples))
							}
					} else {
						if(is.null(columns)) {
							out <- NULL
							return(NULL)
						} else {
							data <- as.data.frame(object@list_of_data[names(object@list_of_data) %in% which_data][[1]])
							phenotype <- object@list_of_row_data[[which_data]]
							out <- phenotype[rownames(phenotype) %in% rownames(data), columns]
							out <- out[order(rownames(out)), ]
							timepoints <- unlist(lapply(strsplit(rownames(out), split="_"), function(x) return(x[2])))
							samples <- unlist(lapply(strsplit(rownames(out), split="_"), function(x) return(x[1])))
							levels <- sort(unique(as.numeric(unlist(lapply(strsplit(timepoints, split="t"), function(x) return(x[2]))))))
							timepoints <- factor(timepoints, levels=paste("t",levels,sep=""))
							out <- as.data.frame(cbind(out, timepoints, samples))
						}
					}
					return(out)
	}) 

#metadata_samples <- get_metadata_for_rows(object=data, which_data="lipid_data", 
#											columns=c("ADNI_MEM", "ADNI_LAN", "ADNI_EF", "APOEGrp", "DXGrp_longi", "PTGENDER", "Age", "BMI"), 
#											names=NULL)

