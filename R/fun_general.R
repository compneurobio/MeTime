# Loading all the packages needed 
#'@import dplyr
#'@import igraph
#'@import visNetwork
#'@import shiny
#'@import shinydashboard
#'@import shinydashboardPlus
#'@import dashboardthemes
#'@import factoextra
#'@import plyr
#'@import DT
#'@import plotly
#'@import M3C
#'@import umap
#'@import ggplot2 
#'@import RColorBrewer
#'@import tidyverse
#'@import GeneNet
#'@import longitudinal
#'@import tidyverse
#'@import glmnet

#' Function to get a palette of distinct colorblind friendly colors.
#' @param n number of colors wanted in the palette
#' @return a color palette vector with colors in the form of hex codes
#' @export 

get_palette <- function(n) {
	#loading the package to get colors 
	#require(RColorBrewer)
	#extracting all the colorblind friendly colors
	colors <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
	col_vec <- unlist(mapply(brewer.pal, colors$maxcolors, rownames(colors)))
	col_vec <- unique(col_vec)
	#get distinct colors by converting the hex to rgb and then to HSL values
	#We then order the colors based on hue and the differences between each are ranked
	#we then retain the 'n' colors needed
	hue_df <- data.frame(color=col_vec, hue=rgb2hsv(col2rgb(col_vec))[1,])#here first row is hue values
	hue_df <- hue_df[order(hue_df$hue),]
  	hue_df <- rbind(hue_df, hue_df[1,])
  	hue_df$hue[nrow(hue_df)] <- hue_df$hue[nrow(hue_df)] + 1
  	hue_df$diff <- c(1, diff(hue_df$hue))
  	hue_df <- hue_df[order(-hue_df$diff),]
  	hue_df <- hue_df[-1, ]
  	return(hue_df$color[1:n])
}

#' Function to split the list of dataframes into a nested list with each dataframe 
#' being split into into dataframes of different timepoints
#' @param object An object of class metab_analyser
#' @return list_of_data with each dataframe being broken into a list of dataframes with respect to the timepoint they belong to
#' @export
setGeneric("split_acc_time", function(object) standardGeneric("split_acc_time") )

setMethod("split_acc_time", "metab_analyser", function(object) {
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
#' @param object An object of class metab_anaylser
#' @param time_splitter A boolean input: True leads to splitting of the data wrt time, 
#'	 					False returns all the dataframes as they are with common rows
#' @return list_of_data with common samples across all time points
#' @export
setGeneric("common_sample_extractor", function(object, time_splitter=FALSE) standardGeneric("common_sample_extractor") )

setMethod("common_sample_extractor", "metab_analyser",function(object, time_splitter=FALSE) {
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

#' Function to Obtain textual information for visualization in interactive plots
#' @param data a dataframe with plotting data along with other variables for visualization
#' @param colnames a character vector with the names of the variables that you want to see on the plot
#' @return a vector with strings that can be parsed into plot_ly text.
#' @export
get_text <- function(data, colnames) {	
		strings_vector <- c()
		count <- 1
		text <- c()
		for(i in 1:length(rownames(data))) {
			for(j in 1:length(colnames)) {
				if(j==1) {
						text <- paste("<br> ", colnames[j], " : ", data[i, colnames[j]], sep="")
					} else {
						text <- paste(text, "<br> ", colnames[j], " : ", data[i, colnames[j]], sep="")
					}
				}
				strings_vector[count] <- text
				count <- count + 1
			} 
		return(strings_vector)
	} 

#' Function to Convert S4 object of class metab_analyser to an S3 object with same architecture
#' @param object An object of class metab_analyser
#' @return An S3 object of the same data as metab_analyser in other words all slots are now converted into nested lists
#' @export
setGeneric("convert_s4_to_s3", function(object) standardGeneric("convert_s4_to_s3"))

setMethod("convert_s4_to_s3", "metab_analyser", function(object) {
		#will add based on analysis - make it module wise or open for suggestions
		return(list(list_of_data=object@list_of_data, list_of_col_data=object@list_of_col_data, list_of_row_data=object@list_of_row_data, annotations=object@annotations))
	})

# #' Function to split data wrt to needed timepoints. Converts S4 back to S3 as well for other analysis.
# #' @description Similar to split_acc_to_time however the timepoints can be manually chosen
# #' @param object An object of class metab_analyser
# #' @param timepoints A vector with timepoints of interest
# #' @return An S3 object(nested list) with the same architecture as that of class metab_analyser
# #' @export
# setGeneric("split_data_wrt_timepoints", function(object, timepoints) standardGeneric("split_data_wrt_timepoints"))
# setMethod("split_data_wrt_timepoints", "metab_analyser", function(object, timepoints) {
# 		data <- split_acc_time(object)
# 		data <- data[names(data) %in% timepoints]
# 		return(list(list_of_data=data,list_of_row_data=object@list_of_row_data, list_of_col_data=object@list_of_col_data))
# 	}) 

#' Function to prepare and preprocess S4 objects to use it for gaussian gaphical models. Also converts S4 to S3
#' @param object An object of class metab_analyser
#' @param which_type two choices either: 1) single -  converts S4 to S3 and returns the nested list
#' 										 2) multi - extracts common samples across the dataframes and returns an S3 nested list
#' @return An S3 object(nested list) with the same architecture as that of class metab_analyser
#' @export 
setGeneric("prep_data_for_ggms", function(object, which_type, mlp_or_temp) standardGeneric("prep_data_for_ggms"))
setMethod("prep_data_for_ggms", "metab_analyser", function(object, which_type, mlp_or_temp) {
		if(which_type %in% "multi") {
			object@list_of_data <- common_sample_extractor(object)
			if(mlp_or_temp) {
				object@list_of_data <- split_acc_time(object)
				object <- convert_s4_to_s3(object)
			} else {
				object <- convert_s4_to_s3(object)
			}
		} else if(which_type %in% "single") {
			if(mlp_or_temp) {
				object@list_of_data <- split_acc_time(object)
				object <- convert_s4_to_s3(object)
			} else {
				object <- convert_s4_to_s3(object)
			}
		} else {
			stop("Check the input for which_type: only allowed inputs and multi and single")
		}
		return(object)
	})

#' Function to add measurements taken at screening time for samples to be added to all timepoints
#' @param object An object of class metab_analyser
#' @param vars A character naming the vars of interest
#' @return phenotype data which can be replaced into the original object or use it separately with a different object
#' @export
setGeneric("add_screening_vars", function(object, vars) standardGeneric("add_screening_vars"))

setMethod("add_screening_vars", "metab_analyser", function(object, vars) {
		phenotype_name <- object@annotations$phenotype
		phenotype <- object@list_of_data[[phenotype_name]]
		screening <- phenotype[grep("-1|-2", rownames(phenotype)), ]
		new_rows <- as.data.frame(screening[, vars])
		new_rows <- na.omit(new_rows)
		sample_names <- unlist(lapply(strsplit(rownames(screening), split="_t"), function(x) return(x[1])))
		for(i in 1:length(sample_names)) {
			index <- grep(sample_names[i], rownames(phenotype)) 
			for(j in index) {
				phenotype[j, vars] <- screening[rownames(screening)[i], vars]
			}
		}
		return(phenotype)
	})

#' Function to know the number of timepoints and the total number of samples available at that point
#' @param object An object of class metab_analyser
#' @param which_data Name of the dataset in context
#' @return A data table with timepoints and number of samples at each timepoint
#' @export
setGeneric("get_samples_and_timepoints", function(object, which_data) standardGeneric("get_samples_and_timepoints"))

setMethod("get_samples_and_timepoints", "metab_analyser", function(object, which_data){
		data <- object@list_of_data[names(object@list_of_data) %in% which_data]
		data <- data[[1]]
		unique_timepoints <- unique(unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[2]))))
		levels <- order(as.numeric(unique(unlist(lapply(strsplit(rownames(data), split="_t"), function(x) return(x[2]))))))
		unique_timepoints <- unique_timepoints[levels]
		samples_count <- vapply(unique_timepoints, function(x) {
						index <- grep(x, rownames(data))
						return(length(index))
					}, numeric(1))
		newdata <- as.data.frame(cbind(as.character(unique_timepoints), samples_count))
		colnames(newdata) <- c("timepoints", "number of samples")
		return(newdata)
	})

#' Function to extract metadata of the metabolites
#' @param object S4 object of class metab_analyse
#' @param which_data choose the dataset from which metabolites will be extracted for metadata
#' @param metab_groups choose the column that has metabolite groups
#' @param metab_ids chodse the column that has metabolite names
#' @return metadata dataframe with names, groups and class
setGeneric("get_metadata_for_plotting", function(object, which_data, metab_groups, metab_ids) standardGeneric("get_metadata_for_plotting"))

setMethod("get_metadata_for_plotting", "metab_analyser", function(object, which_data, metab_groups, metab_ids) {
			col_data_names <- gsub("_data", "_col_data", which_data)
    		row_data_names <- gsub("_data", "_row_data", which_data)
    		list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
    		list_of_col_data <- object@list_of_col_data[names(object@list_of_col_data) %in% col_data_names]
    		list_of_row_data <- object@list_of_row_data[names(object@list_of_row_data) %in% row_data_names]
			if(length(which_data) > 1) {
				list_of_metadata <- list()
				for(i in 1:length(metab_ids)) {
					list_of_metadata[[i]] <- list_of_col_data[[i]][,c(metab_ids[i], metab_groups[i])]
					class <- rep(which_data[i], each=length(list_of_metadata[[i]][,1]))
					list_of_metadata[[i]][,3] <- class
					colnames(list_of_metadata[[i]]) <- c("name", "group", "class")
				}
				metadata <- lapply(list_of_metadata, as.data.frame)
				metadata <- do.call(rbind, metadata)
			} else {
				col_data <- list_of_col_data[[1]]
				class <- rep(which_data, each=col_data[,1])
				metadata <- col_data[ ,c(metab_ids, metab_groups)]
				metadata <- cbind(metadata, class)
				colnames(metadata) <- c("name", "group", "class")
			}
			return(metadata)
	})  

#lol <- get_metadata_for_plotting(object=object, which_data=c("lipid_data", "nmr_data"), metab_ids=c("metabolite", "id"), metab_groups=c("sub_pathway", "Group"))

#' Function to add indices i.e. RID and timepoints to data matrix
#' @description add ids to check for data with ease
#' @param data data matrix with rownames as adni_ids
#' @return data matrix with new added columns
#' @export
adni_add_index <- function(data){
  out = data %>% 
    dplyr::mutate(adni_id = .[] %>% rownames(),
                  rid = .[] %>% rownames() %>% adni_read_id(to="RID") %>% as.numeric(),
                  timepoint = .[] %>% rownames() %>% adni_read_id(to="timepoint") %>% as.numeric())
  return(out)
}

#' Function to remove indices i.e. RID and timepoints from data matrix
#' @description remove ids that are used to check for data with ease(see adni_add_index())
#' @param data data matrix with rownames as adni_ids
#' @return data matrix without index columns
#' @export
adni_rm_index <- function(data){
  rm_col = intersect(names(data), c("adni_id","RID","rid","timepoint","tp","subject"))
  if(length(rm_col)==0){
    out=data
  }else{
    out = data %>% 
      select(-c(rm_col))
    
  }
  return(out)
}

#' Function to extract timepoints of ineterest for longitudnal analysis
#' @description function to extract the timepoints that are of importance in longitudnal analysis similar to split_acc_to_time except is applied on data matrix rather than the object itself
#' @param x.data data matrix of interest
#' @param tp timepoints of interest
#' @return data matrix with only timepoints of interest
#' @export
adni_filter_full_tp<-function(x.data, tp, by=c("timepoint")){
  if(!all(by %in% names(x.data))) stop("adni_tp_overlap() could not match by")
  my_rid <- x.data %>%     
    dplyr::select(all_of(c("RID",by))) %>% 
    dplyr::filter(get(by) %in% tp)%>% 
    dplyr::group_by(RID) %>% 
    dplyr::count(RID) %>% 
    dplyr::filter(n == length(tp))
  
  out= x.data %>% 
    dplyr::filter(timepoint %in% tp,
                  RID %in% my_rid[["RID"]])
  return(out)
}
