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

#' Function to know the number of timepoints and the total number of samples available at that point
#' @description A method applied onto s4 object of class "metab_analyser" so as to obtain the number of unique samples available
#' at each timepoint. 
#' @examples
#' # newdata <- get_samples_and_timepoints(object=metab_analyser_object, which_data="Name of dataset of interest")
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

#' Function to extract metadata of the metabolites and samples for visualization(both dimensionality reduction and GGMs)
#' @description A method applied onto s4 object of class "metab_analyser" so as to obtain metadata of the metabolites and samples.
#' Metadata includes their ontology that is the pathway they belong to and also the class or the dataset type. can also add colors
#' for the metabolites for visualization as a separate column. For samples the metadata is basically the columns of interest from the phenotype table 
#' that can be used to see sample information in the interactive plot.
#' @examples
#' # metadata_list <- get_metadata_for_plotting(object=metab_analyser_object, which_data="name/s of datasets", 
#'											metab_groups="colname/s of the group column in each dataset in order"
#' 											metab_ids="colname/s of the metabolite_name column in each dataset in order", 
#'											cols_for_vis_samples="colnames of phenotype data for samples", screeing_vars=TRUE/FALSE)
#' @param object S4 object of class metab_analyse
#' @param which_data choose the dataset from which metabolites will be extracted for metadata
#' @param metab_groups choose the column that has metabolite groups
#' @param metab_ids chodse the column that has metabolite names
#' @param cols_for_vis_samples character vector representing the name of the columns in phenotype data
#' @param screening_vars character vector representing the measurements obtained at baseline to be added to all timepoints for visualization.
#' is set to NULL if nothing is added to it
#' @return metadata dataframe with names, groups and class
setGeneric("get_metadata_for_plotting", function(object, which_data, metab_groups, metab_ids, cols_for_vis_samples, screening_vars) standardGeneric("get_metadata_for_plotting"))

setMethod("get_metadata_for_plotting", "metab_analyser", function(object, which_data, metab_groups, metab_ids, cols_for_vis_samples, screening_vars=NULL) {
    	#Extracting metadata for metabolites and samples simultaneously
    	list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
    	list_of_col_data <- object@list_of_col_data[names(object@list_of_col_data) %in% which_data]
    	if(!is.null(screening_vars)) {
					object@list_of_data[[object@annotations$phenotype]] <- add_screening_vars(object, vars=screening_vars)
			}
			phenotype <- object@list_of_data[[object@annotations$phenotype]]
			if(length(which_data) > 1) {
				list_of_metadata_metabs <- list()
				for(i in 1:length(metab_ids)) {
					list_of_metadata_metabs[[i]] <- list_of_col_data[[i]][,c(metab_ids[i], metab_groups[i])]
					class <- rep(which_data[i], each=length(list_of_metadata_metabs[[i]][,1]))
					list_of_metadata_metabs[[i]][,3] <- class
					colnames(list_of_metadata_metabs[[i]]) <- c("name", "group", "class")
				}
				metadata_metabs <- lapply(list_of_metadata_metabs, as.data.frame)
				metadata_metabs <- do.call(rbind, metadata_metabs)
				metadata_metabs <- as.data.frame(metadata_metabs)
				metadata_metabs <- metadata_metabs[sort(metadata_metabs$name), ]
				list_of_data <- mod_common_sample_extractor(object)
				list_of_data <- lapply(list_of_data, function(x) return(x[sort(rownames(x)), ]))
				metadata_samples <- as.data.frame(list_of_data[[object@annotations$phenotype]][,cols_for_vis_samples])
			} else {
				col_data <- list_of_col_data[[1]]
				class <- rep(which_data, each=col_data[,1])
				metadata_metabs <- col_data[ ,c(metab_ids, metab_groups)]
				metadata_metabs <- cbind(metadata_metabs, class)
				colnames(metadata_metabs) <- c("name", "group", "class")
				metadata_metabs <- as.data.frame(metadata_metabs)
				metadata_metabs <- metadata_metabs[sort(metadata_metabs$name), ]
				data <- as.data.frame(object@list_of_data[names(object@list_of_data) %in% which_data][[1]])
				phenotype <- object@list_of_data[[object@annotations$phenotype]]
				metadata_samples <- phenotype[rownames(phenotype) %in% rownames(data), cols_for_vis_samples] 
			}
			return(list(metab=metadata_metabs, samples=metadata_samples))
	})  


#' Function to pack all the data into a single object of class "metab_analyser" 
#'
#' @description This function creates an object for MetabAnalyze from a dataset.
#' @examples
#' # new_metab_analyser_object <- get_make_metab_object(data=data_frame, col_data=col_data_frame, row_data=row_data, name="name of the new dataset", 
#'                                annotations_index=list(phenotype="name of phenotype", medication="name of medication"))
#' @param data data.frame containing data 
#' @param col_data data.frame containing col_data: id column of col data has to match colnames of data
#' @param row_data data.frame containing row_data: id column of row data has to match rownames of data
#' @param annotations_index a list to be filled as follows = list(phenotype="Name or index of the file/list", medication="Name or index of the files/list")
#' @param name character. Name you want to assign to the new dataset that is being added on
#' @return An object of class metab_analyser
#' @export

get_make_metab_object <- function(data, col_data, row_data, annotations_index, name=NULL) {
  if(is.null(name)) name <- "set1"
  if(!all(rownames(data) %in% row_data$id) & !all(colnames(data) %in% col_data$id)) stop("id of col or row data do not match dataframe")
  if(!all(c("id","rid","timepoint") %in% names(row_data))) stop("id, subject or timepoint column missing")
  
  list_of_data <- list()
  list_of_data[[name]] <- data
  
  list_of_col_data <- list()
  list_of_col_data[[name]] <- col_data 
  
  list_of_row_data <- list()
  list_of_row_data[[name]] <- row_data
  if(!is.null(annotations_index)) {
  	metab_object <- new("metab_analyser", 
                      list_of_data=list_of_data, 
                      list_of_col_data=list_of_col_data, 
                      list_of_row_data=list_of_row_data,
                      annotations=annotations_index)
  } else {
  	metab_object <- new("metab_analyser", 
                      list_of_data=list_of_data, 
                      list_of_col_data=list_of_col_data, 
                      list_of_row_data=list_of_row_data)
  }
  
  return(metab_object)
}

#' This function appends an object for MetabAnalyze with a new dataset.
#' @description function to apply on metab_analyse object to append a new dataset into the existing object
#' @examples # append data frames into the metab_analyser object
#' appended_object <- get_append_metab_object(object=metab_analyser_object, data=data, row_data=data, col_data=col_data, name="name of the new dataset")
#' @param object S4 MetabAnalyze object
#' @param data data.frame containing data 
#' @param col_data data.frame containing col_data: id column of col data has to match colnames of data
#' @param row_data data.frame containing row_data: id column of row data has to match rownames of data
#' @param name Name of the new dataset
#' @return An object of class metab_analyser
#' @export
setGeneric("get_append_metab_object", function(object, data, col_data, row_data, name) standardGeneric("get_append_metab_object"))
setMethod("get_append_metab_object", "metab_analyser",function(object, data, col_data, row_data, name=NULL) {
  if(is.null(name)) name <- "set1"
  if(!all(rownames(data) %in% row_data$id) & !all(colnames(data) %in% col_data$id)) stop("id of col or row data do not match dataframe")
  if(!all(c("id","rid","timepoint") %in% names(row_data))) stop("id, subject or timepoint column missing")
  
  object@list_of_data[[name]] <- data
  object@list_of_col_data[[name]] <- col_data
  object@list_of_row_data[[name]] <- row_data
  

  return(object)
})

#' Function to pack all the data into a single object of class "metab_analyser" 
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
#' @return An object of class metab_analyser
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
	metab_object <- new("metab_analyser", list_of_data=list_of_data, list_of_col_data=list_of_col_data, 
										list_of_row_data=list_of_row_data,
										annotations=annotations_index)
	return(metab_object)
}

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



#' Function to make a plottable object for viz functions
#' @description function to generate metab_plotter object from plot data and metadata
#' @param data_list list of plotable data
#' @param metadata_list list of metadata for each plot in data list
#' @param plot_type type of the plot you want to build. eg: "box", "dot" etc
#' @param aesthetics aesthetics for the plot object
get_make_plotter_object <- function(data_list, metadata_list, plot_type, aesthetics) {
			data_list <- lapply(data_list, function(x) return(x[sort(rownames(x)), ]))
			plot_data <- list()
			empty_plots <- list()
			count <- 1
			for(i in 1:length(data_list)) {
					plot_data[[count]] <- cbind(data_list[[i]], metadata_list[[i]])
					empty_plots[[count]] <- ggplot(plot_data[[count]])
					count <- count + 1
			}
			if(plot_type %in% "dot") {
					empty_plots <- lapply(empty_plots, function(x) return(x + geom_point()))
			} else if(plot_type %in% "heatmaps") {
					empty_plots <- lapply(empty_plots, function(x) return(x + geom_tile()))
			}
			object <- new("metab_plotter", plot_data=plot_data, plot_parameters=empty_plots)
			return(object)
}
#'add aesthetics to plot and so on in this style 
#'lol <- aes(x=x, y=y)
#'empty_plot$mappings <- lol

#' creating metab_plotter class that converts calculations and metadata as a plotable object to parse into viz_dot_plotter, viz_heatmap_plotter etc
#' Contains slots - plot_data: Dataframe with plotting data and metadata for visualization
#' 				  - plot_parameters: ggplot() object with predefined aesthetics 
#'                - aesthetics: list to define aesthetics. Eg: aesthetics=list(x="colname.x", y="colname.y", color="color", shape="shape")
#'                - the example above will be predefined in all the methods that creates this object. 
#' @rdname metab_plotter
#' @export
setClass("metab_plotter", slots=list(plot_data="list", plot_parameters="list", aesthetics="list"))
