#creating reference metab-analyser class that creates an object with full data

#' Constructor to generate an object of class metab_analyser. 
#' contains slots - list_of_data: For the list of all data matrices.
#'				  - list_of_col_data: list of all the col data files in the same order.
#' 				  - list_of_row_data: list of all the row data files in the same order.
#' 				  - annotations: list with phenotype and medication. Each of which is character that represents 
#'									the name of the aforementioned dataset types.  	
#' 	
#' @rdname metab_analyser
#' @export 
setClass("metab_analyser", slots=list(list_of_data="list", list_of_col_data="list", list_of_row_data="list", 
								 annotations="list")) 


#' Function to split the list of dataframes into a nested list with each dataframe 
#' being split into into dataframes of different timepoints
#' @param object An object of class metab_analyser
#' @return list_of_data with each dataframe being broken into a list of dataframes with respect to the timepoint they belong to
#' @export
setGeneric("mod_split_acc_to_time", function(object) standardGeneric("mod_split_acc_to_time") )

setMethod("mod_split_acc_to_time", "metab_analyser", function(object) {
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
	return(list_of_data_temporals)
})

#' Function to get only common samples from the dataframes in list_of_data 
#' @description A method applied on object of class metab_analyse to extract common samples across datasets. Also has an option to split the data according 
#' timepoints(modify_split_acc_time()).
#' @param object An object of class metab_anaylser
#' @param time_splitter A boolean input: True leads to splitting of the data wrt time, 
#'	 					False returns all the dataframes as they are with common rows
#' @return list_of_data with common samples across all time points
#' @export
setGeneric("mod_common_sample_extractor", function(object, time_splitter=FALSE) standardGeneric("mod_common_sample_extractor") )

setMethod("mod_common_sample_extractor", "metab_analyser",function(object, time_splitter=FALSE) {
		list_of_data <- object@list_of_data
		list_of_names <- lapply(list_of_data, function(x) {
					return(rownames(x))
			})
		common_samples <- Reduce(intersect, list_of_names)
		list_of_data <- lapply(list_of_data, function(x) {
					x <- x[rownames(x) %in% common_samples, ]
					return(x)
			})
		if(time_splitter) {
				split_acc_time(object)
		}
		return(list_of_data)
})

#' Function to Convert S4 object of class metab_analyser to an S3 object with same architecture
#' @description converter function to be applied onto metab_analyse object to convert into a standard list of S3 type.
#' @param object An object of class metab_analyser
#' @return An S3 object of the same data as metab_analyser in other words all slots are now converted into nested lists
#' @export
setGeneric("mod_convert_s4_to_s3", function(object) standardGeneric("mod_convert_s4_to_s3"))

setMethod("mod_convert_s4_to_s3", "metab_analyser", function(object) {
		#will add based on analysis - make it module wise or open for suggestions
		return(list(list_of_data=object@list_of_data, list_of_col_data=object@list_of_col_data, list_of_row_data=object@list_of_row_data, annotations=object@annotations))
	})


#' Function to prepare and preprocess S4 objects to use it for gaussian gaphical models. Also converts S4 to S3
#' @description function to be applied onto metab_analyse object to convert into a standard list of S3 type based on the type of GGM analysis to be performed.
#' 
#' @param object An object of class metab_analyser
#' @param which_type two choices either: 1) single -  converts S4 to S3 and returns the nested list
#' 										 2) multi - extracts common samples across the dataframes and returns an S3 nested list
#' @return An S3 object(nested list) with the same architecture as that of class metab_analyser
#' @export 
setGeneric("mod_prep_data_for_ggms", function(object, which_type, mlp_or_temp) standardGeneric("mod_prep_data_for_ggms"))
setMethod("mod_prep_data_for_ggms", "metab_analyser", function(object, which_type, mlp_or_temp) {
		if(which_type %in% "multi") {
			object@list_of_data <- mod_common_sample_extractor(object)
			if(mlp_or_temp) {
				object@list_of_data <- mod_split_acc_to_time(object)
				object <- mod_convert_s4_to_s3(object)
			} else {
				object <- mod_convert_s4_to_s3(object)
			}
		} else if(which_type %in% "single") {
			if(mlp_or_temp) {
				object@list_of_data <- mod_split_acc_time(object)
				object <- mod_convert_s4_to_s3(object)
			} else {
				object <- mod_convert_s4_to_s3(object)
			}
		} else {
			stop("Check the input for which_type: only allowed inputs and multi and single")
		}
		return(object)
	})

#' Function to apply log transformation
#' @description Function to log transform data
#' @param object An object of class metab_analyser
#' @param which_data Name of the dataset to be used
#' @param base base of log to be used
#' @return An object of class metab_analyser with processed data
#' @export
setGeneric("mod_logtrans", function(object, which_data, base) standardGeneric("mod_logtrans"))
setMethod("mod_logtrans", "metab_analyser",function(object, which_data, base=2) {
  #define data to be processed
  data_position <- which(names(object@list_of_data) %in% which_data)
  object@list_of_data[data_position] = lapply(object@list_of_data[data_position] , log, base=base)
  return(object)
})

#' Function to scale the data
#' @description Functions for scaling
#' @param object An object of class metab_analyser
#' @param which_data Name of the dataset to be used
#' @return An object of class metab_analyser with processed data
#' @export
setGeneric("mod_zscore", function(object, which_data) standardGeneric("mod_zscore"))
setMethod("mod_zscore", "metab_analyser", function(object, which_data) {
  #define data to be processed
  data_position <- which(names(object@list_of_data) %in% which_data)
  data_rownames <- rownames(object@list_of_data[data_position])
  object@list_of_data[data_position] = lapply(object@list_of_data[data_position] , scale, center=TRUE, scale=TRUE)
  object@list_of_data[data_position] = lapply(object@list_of_data[data_position] , as.data.frame)
  rownames(object@list_of_data[data_position]) = data_rownames
  return(object)
})


#' Functions for selecting time points
#' @param object An object of class metab_analyser
#' @param timepoints time points to be selected
#' @param which_data Name of the dataset to be used
#' @param full if TRUE subjects are only selected if measured in all selected time points
#' @return An object of class metab_analyser with processed data
#' @export 

setGeneric("mod_filter_tp", function(object, timepoints, full, which_data) standardGeneric("mod_filter_tp"))
setMethod("mod_filter_tp", "metab_analyser", function(object, timepoints, full=TRUE, which_data) {
  # define data to be processed
  data_position <- which(names(object@list_of_data) %in% which_data)
  
  for(i in data_position){
    keep_id <- object@list_of_row_data[[i]] %>% 
      dplyr::select(id, timepoint, rid) %>% 
      dplyr::mutate(timepoint = as.numeric(timepoint)) %>% 
      dplyr::filter(timepoint %in% timepoints)
    if(full){
      full_rid <- keep_id %>% 
        dplyr::count(rid) %>% 
        dplyr::filter(n==length(timepoints))
      keep_id <- keep_id %>% 
        dplyr::filter(rid %in% full_rid$rid)
      object@list_of_row_data[[i]] = object@list_of_row_data[[i]] %>% 
        dplyr::filter(id %in% keep_id$id)
      object@list_of_data[[i]] = object@list_of_data[[i]][keep_id$id,]
    }else{
      object@list_of_row_data[[i]] = object@list_of_row_data[[i]] %>% 
        dplyr::filter(id %in% keep_id$id)
      object@list_of_data[[i]] = object@list_of_data[[i]][keep_id$id,]
    }
  }
  return(object)
})
