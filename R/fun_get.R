#' Get a palette of "n" distinct colorblind friendly colors 
#' @description Function to get a palette of distinct colorblind friendly colors, the distinctiveness is determined by the difference in their hue values.
#' @examples
#' # colors=get_palette(n=10)
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

#' Function to Obtain textual information for visualization in interactive plots
#' @description a standard function to be applied on data matrices or dataframes with the colnames of interest such that the information from
#' columns is visualized in the interactive plot
#' @examples
#' # text = get_text(data=data.frame, colnames=c("names","of","columns", "of", "interest"))
#' @param data a dataframe with plotting data along with other variables for visualization
#' @param colnames a character vector with the names of the variables that you want to see on the plot
#' @return a vector with strings that can be parsed into plot_ly text.
#' @export
get_text_for_plot <- function(data, colnames) {	
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

#' Function to know the number of timepoints and the total number of samples available at that point
#' @description A method applied onto s4 object of class "metime_analyser" so as to obtain the number of unique samples available
#' at each timepoint. 
#' @examples
#' # newdata <- get_samples_and_timepoints(object=metime_analyser_object, which_data="Name of dataset of interest")
#' @param object An object of class metime_analyser
#' @param which_data Name of the dataset in context
#' @return A data table with timepoints and number of samples at each timepoint
#' @export
setGeneric("get_samples_and_timepoints", function(object, which_data) standardGeneric("get_samples_and_timepoints"))

setMethod("get_samples_and_timepoints", "metime_analyser", function(object, which_data){
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


#' Function to pack all the data into a single object of class "metime_analyser" 
#'
#' @description This function creates an object of class metime_analyser from a dataset.
#' @examples
#' # new_metime_analyser_object <- get_make_metab_object(data=data_frame, col_data=col_data_frame, row_data=row_data, name="name of the new dataset", 
#'                                annotations_index=list(phenotype="name of phenotype", medication="name of medication"))
#' @param data data.frame containing data 
#' @param col_data data.frame containing col_data: id column of col data has to match colnames of data
#' @param row_data data.frame containing row_data: id column of row data has to match rownames of data
#' @param annotations_index a list to be filled as follows = list(phenotype="Name or index of the file/list", medication="Name or index of the files/list")
#' @param name character. Name you want to assign to the new dataset that is being added on
#' @return An object of class metime_analyser
#' @export

get_make_analyser_object <- function(data, col_data, row_data, annotations_index, name=NULL) {
  if(is.null(name)) name <- "set1"
  if(!all(rownames(data) %in% row_data$id) & !all(colnames(data) %in% col_data$id)) stop("id of col or row data do not match dataframe")
  if(!all(c("id","subject","timepoint") %in% names(row_data))) stop("id, subject or timepoint column missing")
  
  list_of_data <- list()
  list_of_data[[name]] <- data
  
  list_of_col_data <- list()
  list_of_col_data[[name]] <- col_data 
  
  list_of_row_data <- list()
  list_of_row_data[[name]] <- row_data
  if(!is.null(annotations_index)) {
  	metab_object <- new("metime_analyser", 
                      list_of_data=list_of_data, 
                      list_of_col_data=list_of_col_data, 
                      list_of_row_data=list_of_row_data,
                      annotations=annotations_index)
  } else {
  	metab_object <- new("metime_analyser", 
                      list_of_data=list_of_data, 
                      list_of_col_data=list_of_col_data, 
                      list_of_row_data=list_of_row_data)
  }
  
  return(metab_object)
}

#' This function appends an object of class metime_analyser with a new dataset.
#' @description function to apply on metime_analyse object to append a new dataset into the existing object
#' @examples # append data frames into the metime_analyser object
#' appended_object <- get_append_metab_object(object=metime_analyser_object, data=data, row_data=data, col_data=col_data, name="name of the new dataset")
#' @param object S4 object of class metime_analyser
#' @param data data.frame containing data 
#' @param col_data data.frame containing col_data: id column of col data has to match colnames of data
#' @param row_data data.frame containing row_data: id column of row data has to match rownames of data
#' @param name Name of the new dataset
#' @return An object of class metime_analyser
#' @export
setGeneric("get_append_analyser_object", function(object, data, col_data, row_data, name) standardGeneric("get_append_analyser_object"))
setMethod("get_append_analyser_object", "metime_analyser",function(object, data, col_data, row_data, name=NULL) {
  if(is.null(name)) name <- "set1"
  if(!all(rownames(data) %in% row_data$id) & !all(colnames(data) %in% col_data$id)) stop("id of col or row data do not match dataframe")
  if(!all(c("id","subject","timepoint") %in% names(row_data))) stop("id, subject or timepoint column missing")
  
  object@list_of_data[[name]] <- data
  object@list_of_col_data[[name]] <- col_data
  object@list_of_row_data[[name]] <- row_data
  

  return(object)
})

#' Function to pack all the data into a single object of class "metime_analyser" 
#'
#' @description This function loads all the files from the parent directory. It assumes a 
#' certain naming pattern as follows: "datatype_None|col|row_data.rds" 
#' Any other naming pattern is not allowed. The function first writes 
#' all files into a list and each type of data is packed into its respective 
#' class i.e. col_data, row_data or data
#'
#' @examples

#' # Input in the parent directory from which the data files are to be extracted along with annotations_index to specify phenotype and medication data

#' get_files_and_names(path=/path/to/parent/directory, annotations_index=list(phenotype="Name of phenotype file", medication="name of phenotype file"))

#' @param path Path to the parent directory
#' @param annotations_index a list to be filled as follows = list(phenotype="Name or index of the files", medication="Name or index of the files")
#' @return An object of class metime_analyser
#' @export

get_files_and_names <- function(path, annotations_index) {
	#path <- input$files$datapath
	path <- list.files(path, pattern="[.rds|.RDS]", full.names=TRUE)	
	data_list <- lapply(path, function(x) {
				x <- readRDS(x)
				return(x)
		})
	col_data_index <- grep("*_col_*", path)
	row_data_index <- grep("*_row_*", path)
	names_list <- lapply(path, function(x) {
					dummy <- unlist(lapply(strsplit(x, split="/"), function(b) return(b[length(b)])))
					dummy <- unlist(lapply(strsplit(dummy, split=".", fixed=TRUE), function(a) return(a[1])))
					return(dummy)
		})	
	list_of_data <- data_list[-c(col_data_index, row_data_index)]
	list_of_col_data <- data_list[col_data_index]
	list_of_row_data <- data_list[row_data_index]
	names(list_of_data) <- names_list[-c(col_data_index, row_data_index)]
	names(list_of_col_data) <- names_list[col_data_index]
	names(list_of_col_data) <- gsub("_col_data", "_data", names(list_of_col_data))
	names(list_of_row_data) <- names_list[row_data_index]
	names(list_of_row_data) <- gsub("_row_data", "_data", names(list_of_row_data))
	metab_object <- new("metime_analyser", list_of_data=list_of_data, list_of_col_data=list_of_col_data, 
										list_of_row_data=list_of_row_data,
										annotations=annotations_index)
	return(metab_object)
}

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



#' Function to make a plottable object for viz functions
#' @description function to generate metime_plotter object from plot data and metadata
#' @param data dataframe of plotable data obtained from calc object
#' @param metadata dataframe with the metadata for the plot table mentioned above. To obtain these see
#' get_metadata_for_rows() and get_metadata_for_columns()
#' @param calc_type A character to specify type of calculation - will be used for comp_ functions
#' @param calc_info A string to define the information about calculation
#' @param plot_type type of the plot you want to build. eg: "box", "dot" etc. Its a character vector
#' @param style Style of plot, accepted inputs are "ggplot", "circos" and "visNetwork". Is a singular option.
#' @export
get_make_plotter_object <- function(data, metadata, calc_type, calc_info, plot_type, style) {
			plot_data <- list()
			empty_plots <- list()
			if(style %in% "visNetwork") {
					plot_data[["node"]] <- data$node
					plot_data[["edge"]] <- data$edge
					plot_data[["metadata"]] <- metadata
			} else {
					data <- data[order(rownames(data)), ]
					metadata <- metadata[rownames(metadata) %in% rownames(data), ]
					plot_data[[1]] <- as.data.frame(cbind(data, metadata))
			}
			for(i in 1:length(plot_type)) {
				if(style %in% "ggplot") {
						empty_plots[[i]] <- ggplot(plot_data[[i]])
				} else if(style %in% "circos") {
						empty_plots[[i]] <- NULL
				} else if(style %in% "visNetwork") {
						empty_plots[[i]] <- NULL
				}
			}
			object <- new("metime_plotter", plot_data=plot_data, plot=empty_plots, calc_type=calc_type, calc_info=calc_info, plot_type=plot_type, style=style)
			return(object)
}


#' creating metime_plotter class that converts calculations and metadata as a plotable object to parse 
#' into viz_plotter
#' Contains slots - plot_data: list of Dataframe(s) with plotting data and metadata for visualization.
#' 									Dataframes is an option only for visNetwork() plots. Need a list of two dataframes: 
#' 									Nodes dataframe and edge dataframe named as .$node and .$edge
#' 				  			- plot: ggplot(), circos() or visNetwork() object  
#'                - calc_type: A vector to specify type of calculation - will be used for comp_ functions
#'                - calc_info: string to define the information about calculation
#' 								- plot_type: A character vector to define the type of plots that are needed.
#' 								- style: Character that defines the style of plot i.e. a ggplot(), circos() or
#'									visNetwork() plot. Is always a singular input. Cannot have two styles in one object.
#' @rdname metime_plotter
#' @export
setClass("metime_plotter", slots=list(plot_data="list", plot="list", calc_type="character", calc_info="character", 
						plot_type="character", style="character"))

##Think of a check function for both classes separately
#duplicated ids - check timpoint and subject 


#' Get metadata for columns(in most cases for metabolites)
#' @description function to generate a metadata list for building the MeTime plotter object
#' @param object S4 object of class MeTime Analyser
#' @param which_data Names of dataset/s to be used
#' @param columns A list of character vectors for the columns of interest. Length of the list should be
#' same as length of which_data
#' @param names A Character vector with the new names for the columns mentioned above
#' @param index_of_names character vector to define the name of the column in which names of the variables are stored
#' @return data.frame with metadata information
#' @export
setGeneric("get_metadata_for_columns", function(object, which_data, columns, names, index_of_names) standardGeneric("get_metadata_for_columns"))
setMethod("get_metadata_for_columns", "metime_analyser", function(object, which_data, columns, names, index_of_names) {
				list_of_col_data <- object@list_of_col_data[names(object@list_of_col_data) %in% which_data]
				if(length(which_data)>1) {
						list_of_metadata_metabs <- list()
						for(i in 1:length(metab_ids)) {
								list_of_metadata_metabs[[i]] <- list_of_col_data[[i]][, columns[[i]]]
								class <- rep(which_data[i], each=length(list_of_metadata_metabs[[i]][,1]))
								list_of_metadata_metabs[[i]] <- list_of_metadata_metabs[[i]][order(list_of_metadata_metabs[[1]][,index_of_names]), ]
								colnames(list_of_metadata_metabs[[i]]) <- names
								list_of_metadata_metabs[[i]] <- as.data.frame(cbind(list_of_metadata_metabs[[i]], class))
						}
						metadata_metabs <- lapply(list_of_metadata_metabs, as.data.frame)
						metadata_metabs <- do.call(rbind, metadata_metabs)
						metadata_metabs <- as.data.frame(metadata_metabs)
						rownames(metadata_metabs) <- metadata_metabs[,index_of_names] 
				} else {
						col_data <- list_of_col_data[[1]]
						class <- rep(which_data, each=length(col_data[,1]))
						metadata_metabs <- col_data[ ,columns[[1]]]
						metadata_metabs <- as.data.frame(cbind(metadata_metabs, class))
						metadata_metabs <- metadata_metabs[order(metadata_metabs[,index_of_names]), ]
						rownames(metadata_metabs) <- metadata_metabs[,index_of_names]
						colnames(metadata_metabs) <- c(names, "class")
				}
				return(metadata_metabs)
	})

#metadata_metabs <- get_metadata_for_columns(object=data, which_data="lipid_data", columns=list(c("metabolite", "sub_pathway")), 
#										names=c("name", "group"), index_of_names="metabolite")


#' Get metadata for rows(in most cases for samples)
#' @description function to generate a metadata list for building the MeTime plotter object
#' @param object S4 object of class MeTime Analyser
#' @param which_data Names of dataset/s to be used
#' @param columns A list of character vectors for the columns of interest. Length of the list should be
#' same as length of which_data
#' @param names A Character vector with the new names for the columns mentioned above
#' @return data.frame with metadata information for rows
#' @export
setGeneric("get_metadata_for_rows", function(object, which_data, columns, names) standardGeneric("get_metadata_for_rows"))
setMethod("get_metadata_for_rows", "metime_analyser", function(object, which_data, columns, names) {
					list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
					if(length(which_data) > 1) {
							list_of_data <- mod_extract_common_samples(object)
							list_of_data <- lapply(list_of_data, function(x) return(x[order(rownames(x)), ]))
							metadata_samples <- as.data.frame(list_of_data[[object@annotations$phenotype]][,columns])
							timepoints <- unlist(lapply(strsplit(rownames(metadata_samples), split="_"), function(x) return(x[2])))
							levels <- sort(unique(as.numeric(unlist(lapply(strsplit(timepoints, split="t"), function(x) return(x[2]))))))
							timepoints <- factor(timepoints, levels=paste("t",levels,sep=""))
							metadata_samples <- as.data.frame(cbind(metadata_samples, timepoints))
					} else {
							data <- as.data.frame(object@list_of_data[names(object@list_of_data) %in% which_data][[1]])
							phenotype <- object@list_of_data[[object@annotations$phenotype]]
							metadata_samples <- phenotype[rownames(phenotype) %in% rownames(data), columns]
							metadata_samples <- metadata_samples[order(rownames(metadata_samples)), ]
							timepoints <- unlist(lapply(strsplit(rownames(metadata_samples), split="_"), function(x) return(x[2])))
							levels <- sort(unique(as.numeric(unlist(lapply(strsplit(timepoints, split="t"), function(x) return(x[2]))))))
							timepoints <- factor(timepoints, levels=paste("t",levels,sep=""))
							metadata_samples <- as.data.frame(cbind(metadata_samples, timepoints))
					}
					return(metadata_samples)
	}) 

#metadata_samples <- get_metadata_for_rows(object=data, which_data="lipid_data", 
#											columns=c("ADNI_MEM", "ADNI_LAN", "ADNI_EF", "APOEGrp", "DXGrp_longi", "PTGENDER", "Age", "BMI"), 
#											names=NULL)