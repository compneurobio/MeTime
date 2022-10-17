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


#' Function to split data acoording to time
#' @description Function to split the list of dataframes into a nested list with each dataframe 
#' being split into into dataframes of different timepoints
#' @examples #splitting data according to time
#' new_data <- mod_split_acc_to_time(object=metime_analyser_object)
#' @param object An object of class metime_analyser
#' @return list_of_data with each dataframe being broken into a list of dataframes with respect to the timepoint they belong to
#' @export
setGeneric("mod_split_acc_to_time", function(object) standardGeneric("mod_split_acc_to_time") )

setMethod("mod_split_acc_to_time", "metime_analyser", function(object) {
	list_of_data <- object@list_of_data
	list_of_data_temporals <- lapply(list_of_data, function(data) {
		names <- strsplit(rownames(data), split="_")
		times <- unlist(lapply(names, function(x) return(x[2])))
		indList <- split(seq_along(times), times)

		#timepoint separation of the data

		list_of_temporals <- lapply(indList, function(x) {
				return(data[x,])
		})
		names(list_of_temporals) <- names(indList)
		return(list_of_temporals)
	})
	names(list_of_data_temporals) <- names(list_of_data)
	object@list_of_data <- list_of_data_temporals
	out <- object
	return(out)
})

#' Function to get only common samples from the dataframes in list_of_data 
#' @description A method applied on object of class metime_analyse to extract common samples across datasets. Also has an option to split the data according 
#' timepoints(mod_split_acc_time()).
#' @examples # extracting common samples across all datasets
#' new_list_of_data <- mod_common_sample_extractor(object=metime_analyser_object)
#' @param object An object of class metime_anaylser
#' @param time_splitter A boolean input: True leads to splitting of the data wrt time, 
#'	 					False returns all the dataframes as they are with common rows
#' @return list_of_data with common samples across all time points
#' @export
setGeneric("mod_extract_common_samples", function(object, time_splitter=FALSE) standardGeneric("mod_extract_common_samples") )

setMethod("mod_extract_common_samples", "metime_analyser",function(object, time_splitter=FALSE) {
		list_of_names <- lapply(object@list_of_data, function(x) {
					return(rownames(x))
			})
		common_samples <- Reduce(intersect, list_of_names)
		object@list_of_data <- lapply(object@list_of_data, function(x) {
					x <- x[rownames(x) %in% common_samples, ]
					x <- x[order(rownames(x)), ]
					return(x)
			})
		if(time_splitter) {
				object <- mod_split_acc_to_time(object)
		}
		out <- object
		return(out)
})

#' Function to Convert S4 object of class metime_analyser to an S3 object with same architecture
#' @description converter function to be applied onto metime_analyse object to convert into a standard list of S3 type.
#' @examples # convert S4 object to a list
#' s3_list <- mod_convert_s4_to_s3(object=metime_analyser_object)
#' @param object An object of class metime_analyser
#' @return An S3 object of the same data as metime_analyser in other words all slots are now converted into nested lists
#' @export
setGeneric("mod_convert_s4_to_s3", function(object) standardGeneric("mod_convert_s4_to_s3"))

setMethod("mod_convert_s4_to_s3", "metime_analyser", function(object) {
		#will add based on analysis - make it module wise or open for suggestions
		out <- list(list_of_data=object@list_of_data, list_of_col_data=object@list_of_col_data, list_of_row_data=object@list_of_row_data, annotations=object@annotations)
		return(out)
	})


#' Function to apply log transformation
#' @description Function to log transform data
#' @examples # example to apply log transformation
#' object <- mod_logtrans(object, which_data="name of the dataset", base=2)
#' @param object An object of class metime_analyser
#' @param which_data Name of the dataset to be used
#' @param base base of log to be used
#' @return An object of class metime_analyser with processed data
#' @export
setGeneric("mod_trans_log", function(object, which_data, base) standardGeneric("mod_trans_log"))
setMethod("mod_trans_log", "metime_analyser",function(object, which_data, base=2) {
  #define data to be processed
  data_position <- which(names(object@list_of_data) %in% which_data)
  object@list_of_data[data_position] = lapply(object@list_of_data[data_position] , log, base=base)
  out <- object
  return(out)
})

#' Function to scale the data
#' @description Functions for scaling
#' @examples # example to apply scaling
#' object <- mod_zscore(object, which_data="name of the dataset")
#' @param object An object of class metime_analyser
#' @param which_data Name of the dataset to be used
#' @return An object of class metime_analyser with processed data
#' @export
setGeneric("mod_trans_zscore", function(object, which_data) standardGeneric("mod_trans_zscore"))
setMethod("mod_trans_zscore", "metime_analyser", function(object, which_data) {
  #define data to be processed
  data_position <- which(names(object@list_of_data) %in% which_data)
  data_rownames <- rownames(object@list_of_data[data_position])
  object@list_of_data[data_position] = lapply(object@list_of_data[data_position] , scale, center=TRUE, scale=TRUE)
  object@list_of_data[data_position] = lapply(object@list_of_data[data_position] , as.data.frame)
  rownames(object@list_of_data[data_position]) = data_rownames
  out <- object
  return(out)
})


#' Functions for selecting time points
#' @description a method applied onto class metime_analyser in order to extract timepoints of interest from a dataset
#' @examples #example to use this function
#' object <- mod_filter_tp(object, timepoints=c("t0","t12","t24"), full=TRUE, which_data="Name of the dataset")
#' @param object An object of class metime_analyser
#' @param timepoints time points to be selected. 
#' @param which_data Name of the dataset to be used
#' @param full if TRUE subjects are only selected if measured in all selected time points
#' @return An object of class metime_analyser with processed data
#' @export 

setGeneric("mod_filter_tp", function(object, timepoints, full, which_data) standardGeneric("mod_filter_tp"))
setMethod("mod_filter_tp", "metime_analyser", function(object, timepoints, full=TRUE, which_data) {
  # define data to be processed
  data_position <- which(names(object@list_of_data) %in% which_data)
  
  for(i in data_position){
    keep_id <- object@list_of_row_data[[i]] %>% 
      dplyr::select(id, time, subject) %>%  
      dplyr::filter(time %in% timepoints)
    if(full) {
      full_rid <- keep_id %>% 
        dplyr::count(subject) %>% 
        dplyr::filter(n==length(timepoints))
      keep_id <- keep_id %>% 
        dplyr::filter(subject %in% full_rid$subject)
      object@list_of_row_data[[i]] = object@list_of_row_data[[i]] %>% 
        dplyr::filter(id %in% keep_id$id)
      object@list_of_data[[i]] = object@list_of_data[[i]][keep_id$id,]
    } else {
      object@list_of_row_data[[i]] = object@list_of_row_data[[i]] %>% 
        dplyr::filter(id %in% keep_id$id)
      object@list_of_data[[i]] = object@list_of_data[[i]][keep_id$id,]
    }
  }
  out <- object
  return(out)
})


#' Function to merge one or more metime_analyser objects
#' @description function to merge multiple metime_analyser objects
#' @param list_of_objects list of metime analyser objects that are to be merged
#' @param annotations_index new list with annotations_index. Can also set to be NULL.
#' @returns A merged metime_analyser object
#' @export
setGeneric("mod_merge_metime_analysers", function(list_of_objects, annotations_index) standardGeneric("mod_merge_metime_analysers"))
setMethod("mod_merge_metime_analysers", "metime_analyser", function(list_of_objects, annotations_index) {
				list_of_data <- lapply(list_of_objects, function(x) return(x@list_of_data))
				list_of_col_data <- lapply(list_of_objects, function(x) return(x@list_of_col_data))
				list_of_row_data <- lapply(list_of_objects, function(x) return(x@list_of_row_data))
				if(is.null(annotations_index)) {
						final_object <- new("metime_analyser", list_of_data=list_of_data, 
							list_of_row_data=list_of_row_data, list_of_col_data=list_of_col_data, annotations=NULL)
				} else {
						final_object <- new("metime_analyser", list_of_data=list_of_data, 
							list_of_row_data=list_of_row_data, list_of_col_data=list_of_col_data, annotations=annotations_index)	
				}
				out <- final_object
				return(out)
	}) 

#' Function to remove NA's from data matrices
#' @description A method applied on S4 object to remove NA's and change data accordingly
#' @param object S4 object of class metime_analyser
#' @param which_data dataset/s for which the method is to be applied
#' @return S4 object with NA's removed and data manipulated accordingly
#' @export
setGeneric("mod_remove_nas", function(object, which_data) standardGeneric("mod_remove_nas"))
setMethod("mod_remove_nas", "metime_analyser", function(object, which_data) {
				list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
				list_of_data <- lapply(list_of_data, function(x) return(na.omit(x)))
				for(i in 1:length(list_of_data)) {
						rows <- object@list_of_row_data[[which_data[i]]]
						cols <- object@list_of_col_data[[which_data[i]]]
						rows <- rows[rownames(rows) %in% rownames(list_of_data[[i]]), ]
						cols <- cols[rownames(cols) %in% colnames(list_of_data[[i]]), ]
						object@list_of_row_data[[i]] <- rows
						object@list_of_col_data[[i]] <- cols					
				}
				out <- object
				return(out)
	})

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

#' Function to remove duplicates
#' @description Function to remove duplicates from the analyser object
#' @param object An S4 object of class metime_analyser
#' @return object after removing duplicated data
#' @export
setGeneric("mod_remove_duplicates", function(object) standardGeneric("mod_remove_duplicates"))
setMethod("mod_remove_duplicates", "metime_analyser", function(object) {
				out <- lapply(names(object@list_of_data), function(x) {
								data <- object@list_of_data[[x]]
								rowdata <- object@list_of_row_data[[x]]
								coldata <- object@list_of_col_data[[x]]
								data <- data[!duplicated(data), ]
								rowdata <- rowdata[!duplicated(rowdata), ]
								coldata <- coldata[!duplicated(coldata), ]
								object@list_of_data[[x]] <- data
								object@list_of_row_data[[x]] <- rowdata
								object@list_of_col_data[[x]] <- coldata
								return(object)
					})
				return(out)
	})



#' Function to stratify data in the metime analyser object
#' @description Function to stratify the data of interest into different objects that can be used
#' to perform calculations according the said stratification variable
#' @param object S4 object of class metime_analyser
#' @param which_data Dataset/datasets to be used for stratification
#' @param variable Phenotype based on which the stratification would be performed
#' @return list of metime_analyser objects which are stratified based on the variable chosen
#' @export
setGeneric("mod_stratify_analyser", function(object, which_data, variable) standardGeneric("mod_stratify_analyser"))
setMethod("mod_stratify_analyser", "metime_analyser", function(object, which_data, variable) {
				out <- lapply(which_data, function(x) {
								data <- object@list_of_data[[x]]
								row_data <- object@list_of_row_data[[x]]
								strat <- as.data.frame(nnet::class.ind(row_data[ ,variable]))
								rownames(strat) <- rownames(row_data)
								strat_list <- lapply(colnames(strat), function(y) {
												strat_samples <- ifelse(strat[,y]==1, rownames(strat), NA)
												strat_samples <- na.omit(strat_samples)
												strat_data <- data[rownames(data) %in% strat_samples, ]
												strat_row_data <- row_data[rownames(row_data) %in% strat_samples, ]
												object@list_of_data[[x]] <- strat_data
												object@list_of_row_data[[x]] <- strat_row_data
												return(object)
									})
								return(strat_list)
					})
				return(out)
	})