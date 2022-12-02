
#' Function for Plotting distributions of phenotypic variables 
#' @description A method to be applied onto s4 object so as to obtain distributions of various phenotypic variables
#' @examples # extracting distribuiton of Age from dataset1
#' plot <- viz_distribution_plotter(object, colname="Age", which_data="dataset1", strats="additional columns for facet wrapping")
#' @param object An object of class metime_analyser
#' @param colname Name of the variable whose distribution is of interest
#' @param which_data Name of the dataset from which the samples will be extracted
#' @param strats Character vector with colnames that are to be used for stratification
#' @param phenotype Logical. If true data will be collected from phenotype_data matrix else from row data 
#' @return a list with either 1) density plot, mean table acc to timepoint and variable type or 
#'								2) bar plot, line plot, and variable type
#' @export
setGeneric("viz_distribution_plotter", function(object, colname, which_data, strats, phenotype) standardGeneric("viz_distribution_plotter") )

setMethod("viz_distribution_plotter", "metime_analyser",function(object, colname, which_data, strats=NULL, phenotype) {
	data <- object@list_of_data[names(object@list_of_data) %in% which_data]
	data <- data[[1]]
	if(phenotype) {
		phenotype <- object@list_of_data[[object@annotations$phenotype]]
	} else {
		phenotype <- object@list_of_row_data[[which_data]]
	}
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
					 geom_line(aes_string(color=colname)) + geom_point(aes_string(color=colname)) + scale_color_manual(values=palette_line) + theme_classic()
		return(out=list(bar_plot=bar_plot, line_plot=line_plot, var_type=var_type))

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
		out <- list(table=mu, density_plot=density_plot, var_type=var_type)
		return(out)
	}
})


