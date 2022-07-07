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


#' Function for Plotting distributions of phenotypic variables 
#' @param object An object of class metab_analyser
#' @param colname Name of the variable whose distribution is of interest
#' @param which_data Name of the dataset from which the samples will be extracted
#'
#' @return a list with either 1) density plot, mean table acc to timepoint and variable type or 
#'								2) bar plot, line plot, and variable type
#' @export
setGeneric("viz_distribution_plotter", function(object, colname, which_data, strats) standardGeneric("viz_distribution_plotter") )

setMethod("viz_distribution_plotter", "metab_analyser",function(object, colname, which_data, strats=NULL) {
	data <- object@list_of_data[names(object@list_of_data) %in% which_data]
	data <- data[[1]]
	phenotype_name <- object@annotations$phenotype
	phenotype <- object@list_of_data[[phenotype_name]]
	var_type <- c()
	vec <- phenotype[,colname]
	vec <- na.omit(vec)
	if(is.character(vec)==TRUE) var_type <- "bar"
	if(is.numeric(vec)==TRUE) {
		if(length(unique(vec)) <= 12) {
				var_type <- "bar"
		} else {
				var_type <- "density"
		}
	}
	if(var_type %in% "bar") {
		phenotype <- phenotype[rownames(phenotype) %in% rownames(data),]
		vec <- phenotype[, colname]
        palette_line <- get_palette(length(unique(vec)))
		palette_timepoints <- get_palette(30)[c(16,2,3,6,7,8,10,13,19,20,23,27,30)]
		timepoints <- unlist(lapply(strsplit(rownames(phenotype), split="_"), function(x) return(x[2])))
		levels <- sort(unique(as.numeric(unlist(lapply(strsplit(timepoints, split="t"), function(x) return(x[2]))))))
		timepoints <- factor(timepoints, levels=paste("t",levels,sep=""))
		plot_data <- data.frame(vec=as.character(vec), timepoints=timepoints)
		plot_data <- table(plot_data)
		plot_data <- reshape2::melt(plot_data)
		colnames(plot_data) <- c(colname, "Timepoints", "Frequency")
		plot_data[,colname] <- as.factor(plot_data[,colname])
		bar_plot <- ggplot(data=plot_data, aes_string(x=colname, y="Frequency", fill="Timepoints")) +
					geom_bar(stat="identity") + scale_fill_manual(values=palette_timepoints)+ theme_classic()
		line_plot <- ggplot(plot_data, aes_string(x="Timepoints", y="Frequency", group=colname)) + 
					 geom_line(aes_string(color=colname)) + geom_point(aes_string(color=colname)) + scale_color_manual(values=palette_line) + theme_minimal()
		return(list_of_plots=list(bar_plot=bar_plot, line_plot=line_plot, var_type=var_type))

	} else {
		phenotype <- phenotype[rownames(phenotype) %in% rownames(data),]
		vec <- phenotype[,colname]
		palette_timepoints <- get_palette(30)[c(16,2,3,6,7,8,10,13,19,20,23,27,30)]
		timepoints <- unlist(lapply(strsplit(rownames(phenotype), split="_"), function(x) return(x[2])))
		if(is.null(strats)) {
			plot_data <- as.data.frame(cbind(as.character(vec), timepoints))
			plot_data[,1] <- as.numeric(plot_data[,1])
			levels <- sort(unique(as.numeric(unlist(lapply(strsplit(timepoints, split="t"), function(x) return(x[2]))))))
			plot_data$timepoints <- factor(timepoints, levels=paste("t", levels, sep=""))
			colnames(plot_data) <- c(colname, "Timepoints")
			plot_data <- na.omit(plot_data)
			mu <- as.data.frame(aggregate(plot_data[,1], list(Timepoints=plot_data[,2]), FUN=mean))
			mu$Timepoints <- factor(mu$Timepoints, levels=paste("t", levels, sep=""))
			colnames(mu)[2] <- "mean"
			density_plot <- ggplot(plot_data, aes_string(x=colname, fill="Timepoints")) + geom_density(alpha=0.3) + 
								geom_vline(data=mu, aes(xintercept=mean, color=Timepoints), linetype="dashed") +
			 					scale_fill_manual(values=palette_timepoints) + theme_classic()
		} else {
			strats_df <- phenotype[,strats]
			plot_data <- as.data.frame(cbind(as.character(vec), timepoints, strats_df))
			plot_data[,1] <- as.numeric(plot_data[,1])
			levels <- sort(unique(as.numeric(unlist(lapply(strsplit(timepoints, split="t"), function(x) return(x[2]))))))
			plot_data$timepoints <- factor(timepoints, levels=paste("t", levels, sep=""))
			colnames(plot_data) <- c(colname, "Timepoints", strats)
			plot_data <- na.omit(plot_data)
			mu <- as.data.frame(aggregate(plot_data[,1], list(Timepoints=plot_data[,2]), FUN=mean))
			mu$Timepoints <- factor(mu$Timepoints, levels=paste("t", levels, sep=""))
			colnames(mu)[2] <- "mean"
			density_plot <- ggplot(plot_data, aes_string(x=colname, fill="Timepoints")) + geom_density(alpha=0.3) + 
								geom_vline(data=mu, aes(xintercept=mean, color=Timepoints), linetype="dashed") +
			 					scale_fill_manual(values=palette_timepoints) + theme_classic() + facet_wrap(strats)
		}
		
		#render table for mean values
		return(list(table=mu, density_plot=density_plot, var_type=var_type))
	}
})

#' Function to dot plot any kind of dot_plotter including for dimensionality reduction
#' @description General function to be implemented on data_list that is obtained after applying a dimensionality reduction method
#' @param data_list list obtained after applying calc_dimensionality_reduction() on metab_analyse object
#' @param metadata_list list obtained after applying get_metadata_for_plotting() on metab_analyse object
#' @param axes_labels character vector to specify the labels of the axes in the order x and y.
#' @param title_metabs character to specify the title of the plot of metabolites
#' @param title_metabs character to specify the title of the plot of metabolites
#' @return a list with both the plots of samples and metabolites. Can be accessed by using ".$samples" and ".$metabs"
viz_dimensionality_reduction <- function(data_list, metadata_list, axes_labels, title_metabs, title_samples) {
	palette <- get_palette(30)
	if(is.null(metadata_list)) {
		return(list(metabs=plot(data_list$metabs), samples=plot(data_list$samples)))
	} else {
		plot <- list()
		data_list <- lapply(data_list, function(x) return(x[sort(rownames(x)),]))
		plot_data_metabs <- as.data.frame(cbind(data_list$metabs, metadata_list$metabs))
		plot_data_samples <- as.data.frame(cbind(data_list$samples, metadata_list$samples))
		plot_metabs <- ggplot(plot_data_metabs, aes_string(x=colnames(plot_data_metabs)[1], y=colnames(plot_data_metabs)[2], 
							color="groups", shape="class", text=paste("metab_name :", plot_data_metabs$name, sep="")))  
						+ geom_point() + labs(x=axes_labels[1], y=axes_labels[2], subtitle=title_metabs) 
		plot_metabs <- ggplotly(metabs)
		timepoints <- unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[2])))
		levels <- unique(sort(as.numeric(unlist(lapply(strsplit(timepoints, split="t", fixed=TRUE), function(x) return(x[2]))))))
		plot_data_samples$timepoints <- factor(timepoints, levels=paste("t", levels, sep=""))
		plot_samples <- ggplot(plot_data_samples, aes_string(x=colnames(plot_data_samples)[1], y=colnames(plot_data_samples)[2], 
						color="timepoints", 
						text=get_text(data=plot_data_samples, colnames=colnames(plot_data_samples)[3:length(colnames(plot_data_samples))]))) +
						geom_point() + labs(x=axes_labels[1], y=axes_labels[2], subtitle=title_samples)
		plot_samples <- ggplotly(plot_samples)
		return(list(metabs=plot_metabs, samples=plot_samples))
	}
}

#viz_conversation_index()