
#Functions for study characteristics


#' Function for Plotting distributions of phenotypic variables 
#' @param object An object of class metab_analyser
#' @param colname Name of the variable whose distribution is of interest
#' @param which_data Name of the dataset from which the samples will be extracted
#'
#' @return a list with either 1) density plot, mean table acc to timepoint and variable type or 
#'								2) bar plot, line plot, and variable type
#' @export
setGeneric("distribution_plotter", function(object, colname, which_data, strats) standardGeneric("distribution_plotter") )

setMethod("distribution_plotter", "metab_analyser",function(object, colname, which_data, strats=NULL) {
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

#sex <- distribution_plotter(object=data, colname="PTGENDER", "")

#' Function to Plot PCA for one dataset with samples as data points
#' @param object An object of class metab_analyser
#' @param cols_for_vis choose the colnames from phenotype data that you want to show in the text labels
#' @param which_data Name of the dataset from which the samples will be extracted
#' 
#' @return an interactive pca plot with text that can be modified.
#' @export
setGeneric("pca_plotter_general", function(object, which_data, cols_for_vis) standardGeneric("pca_plotter_general") )

setMethod("pca_plotter_general", "metab_analyser", function(object, which_data, cols_for_vis) {
	data <- object@list_of_data
	data <- data[names(data) %in% which_data]
	data <- as.data.frame(data[[1]])
	data <- data[, which(apply(data,2,var) !=0)]
	data <- data[which(apply(data,1,var) != 0), ]
	data <- na.omit(data)
	phenotype_name <- object@annotations$phenotype
	phenotype <- object@list_of_data[[phenotype_name]]
	phenotype <- phenotype[rownames(phenotype) %in% rownames(data),]
	timepoints <- unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[2])))
	levels <- unique(sort(as.numeric(unlist(lapply(strsplit(timepoints, split="t", fixed=TRUE), function(x) return(x[2]))))))
	timepoints <- factor(timepoints, levels=paste("t", levels, sep=""))
	pca_individuals <- prcomp(data, scale.=T, center=T)
	#pca_metabolites <- prcomp(t(data), scale.=T, center=T)
	#require(plotly)
	#list_of_pcs <- list(individuals=pca_individuals, metabotlites=metabolites)

	pca_individuals <- as.data.frame(pca_individuals$x[,1:2])
	pca_individuals <- pca_individuals[order(match(rownames(pca_individuals), rownames(phenotype))), , drop = FALSE]
	plot_data <- as.data.frame(cbind(pca_individuals, timepoints, phenotype[,cols_for_vis]))
	colnames(plot_data)[1:2] <- c("PC1", "PC2")
	
	plot_individuals <- plot_ly(plot_data, x=~PC1,y=~PC2, color=~timepoints, colors=get_palette(30)[c(16,2,3,6,7,8,10,13,19,20,23,27,30)], type = 'scatter', 
															mode = 'markers', hoverinfo = 'text', 
														text = get_text(data=plot_data, colnames=cols_for_vis))
	
	return(plot_individuals)
})

#' Function to plot tsne plots for one data set
#' @param object An object of class metab_analyser
#' @param cols_for_vis choose the colnames from phenotype data that you want to show in the text labels
#' @param which_data Name of the dataset from which the samples will be extracted
#' @param metab_ids name of the column with metabolite names in the col data matix
#' @param metab_groups name of the column which has the pathway information of the metabolites
#' 
#' @return a list with two plot objects 1) samples - tSNE plot of the individuals(samples)
#'										 2) metabs - tSNE plot of the metabolites(metabs)
#' @export
setGeneric("tsne_plotter_general", function(object, which_data, metab_ids, metab_groups, cols_for_vis) standardGeneric("tsne_plotter_general"))

setMethod("tsne_plotter_general", "metab_analyser", function(object, which_data, metab_ids, metab_groups, cols_for_vis) {
	data <- object@list_of_data[names(object@list_of_data) %in% which_data]
	data <- data[[1]]
	phenotype_name <- object@annotations$phenotype
	phenotype <- object@list_of_data[[phenotype_name]]
	data <- as.data.frame(data)
	col_data <- object@list_of_col_data[names(object@list_of_col_data) %in% which_data]
	col_data <- as.data.frame(col_data[[1]])
	groups_metabolites <- as.vector(col_data[ ,metab_groups])
	data <- na.omit(data)
	col_data <- col_data[col_data[,metab_ids] %in% colnames(data), ]
	data <- data[ ,match(colnames(data), col_data[ ,metab_ids])]
	phenotype <- phenotype[rownames(phenotype) %in% rownames(data),]
	phenotype <- phenotype[match(rownames(phenotype), rownames(data)),]	
	tsne_plot_metabs <- tsne(data, labels=as.factor(groups_metabolites))
	tsne_samples <- tsne(t(data))
	timepoints <- unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[2])))
	levels <- unique(sort(as.numeric(unlist(lapply(strsplit(timepoints, split="t", fixed=TRUE), function(x) return(x[2]))))))
	tsne_data <- as.data.frame(cbind(as.data.frame(tsne_samples$data), phenotype[,cols_for_vis]))
	tsne_data$timepoints <- factor(timepoints, levels=paste("t", levels, sep=""))
	tsne_plot_samples <- plot_ly(tsne_data, x=~X1, y=~X2, color=~timepoints, colors=get_palette(30)[c(16,2,3,6,7,8,10,13,19,20,23,27,30)], type = 'scatter', 
															mode = 'markers', hoverinfo = 'text', 
														text = get_text(data=tsne_data, colnames=cols_for_vis))
	return(list(metabs=tsne_plot_metabs, samples=tsne_plot_samples))
})

#lol <- tsne_plotter_general(object=data, which_data="lipid_data", metab_groups="sub_pathway", metab_ids="metabolite", cols_for_vis=cols_for_pca_plot)

#UMAP general plotter similar idea as that of tsne_general_plotter

#' Function to plot UMAP plots for one data set
#' @param object An object of class metab_analyser
#' @param cols_for_vis choose the colnames from phenotype data that you want to show in the text labels
#' @param which_data Name of the dataset from which the samples will be extracted
#' @param metab_ids name of the column with metabolite names in the col data matix
#' @param metab_groups name of the column which has the pathway information of the metabolites
#' @return a list with two plot objects 1) samples - UMAP plot of the individuals(samples)
#'										 2) metabs - UMAP plot of the metabolites(metabs)
#' @export

setGeneric("umap_plotter_general", function(object, which_data, metab_ids, metab_groups, cols_for_vis, phenotype_index) standardGeneric("umap_plotter_general")) 

setMethod("umap_plotter_general", "metab_analyser", function(object, which_data, metab_ids, metab_groups, cols_for_vis, phenotype_index) {
	data <- object@list_of_data[names(object@list_of_data) %in% which_data]
	data <- data[[1]]
	phenotype_name <- object@annotations$phenotype
	phenotype <- object@list_of_data[[phenotype_name]]
	data <- as.data.frame(data)
	col_data <- object@list_of_col_data[names(object@list_of_col_data) %in% which_data]
	col_data <- as.data.frame(col_data[[1]])
	groups_metabolites <- as.vector(col_data[ ,metab_groups])
	data <- na.omit(data)
	col_data <- col_data[col_data[,metab_ids] %in% colnames(data),]
	data <- data[ ,match(colnames(data), col_data[ ,metab_ids])]
	phenotype <- phenotype[rownames(phenotype) %in% rownames(data),]
	phenotype <- phenotype[match(rownames(phenotype), rownames(data)),]
	umap_fit <- umap::umap(t(data))
	plot_data <- cbind(as.data.frame(umap_fit$layout), groups_metabolites)
	colnames(plot_data) <- c("UMAP1", "UMAP2", "groups")
	umap_plot_metabs <- ggplot(plot_data, aes(x=UMAP1, y=UMAP2, color=groups)) + geom_point() + labs(x="UMAP1", y="UMAP2", subtitle="UMAP plot for metabolites")
	umap_fit_samples <- umap::umap(data)
	timepoints <- unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[2])))
	levels <- unique(sort(as.numeric(unlist(lapply(strsplit(timepoints, split="t", fixed=TRUE), function(x) return(x[2]))))))
	timepoints <- factor(timepoints, levels=paste("t", levels, sep=""))
	plot_data <- as.data.frame(cbind(as.data.frame(umap_fit_samples$layout), phenotype[,cols_for_vis]))
	plot_data$timepoints <- factor(timepoints, levels=paste("t", levels, sep=""))
	colnames(plot_data)[1:2] <- c("UMAP1","UMAP2")
	umap_plot_samples <- plot_ly(plot_data, x=~UMAP1, y=~UMAP2, color=~timepoints, colors=get_palette(30)[c(16,2,3,6,7,8,10,13,19,20,23,27,30)], type = 'scatter', 
															mode = 'markers', hoverinfo = 'text', 
														text = get_text(data=plot_data, colnames=cols_for_vis))
	return(list(metabs=umap_plot_metabs, samples=umap_plot_samples))
})

#lul <- umap_plotter_general(object=data, which_data="lipid_data", metab_groups="sub_pathway", metab_ids="metabolite", cols_for_vis=cols_for_pca_plot)

#' Function to plot PCA of the individuals and metabolites of common data
#' @param object An object of class metab_analyser
#' @param cols_for_vis choose the colnames from phenotype data that you want to show in the text labels
#' @param which_data a character vector - Names of the dataset from which the samples will be extracted
#' @param metab_ids a character vector - names of the column with metabolite names in the col data matix(please make sure they are in the same order as the above ones)
#' @param metab_groups a character vector - names of the column which has the pathway information of the metabolites(please make sure they are in the same order as the above ones)
#' @param phenotype_index index of the phenotype data. Can input either the name of the phenotype dataset or the index of the same
#' @return a list with two plot objects 1) samples - PCA plot of the individuals(".$samples")
#'										 2) metabs - PCA plot of the metabolites(".$metabs")
#' @export
setGeneric("pca_plotter_wrt_common", function(object, which_data, metab_groups, metab_ids, cols_for_vis, phenotype_index) standardGeneric("pca_plotter_wrt_common"))

setMethod("pca_plotter_wrt_common", "metab_analyser", function(object, which_data, metab_groups, metab_ids, cols_for_vis, phenotype_index) {
		new_list_of_data <- common_sample_extractor(object)
		new_list_of_data <- new_list_of_data[names(new_list_of_data) %in% which_data]
		new_data_for_pca <- do.call(cbind, new_list_of_data)
		new_data_for_pca <- na.omit(new_data_for_pca)
		new_data_for_pca <- as.data.frame(new_data_for_pca)
		data <- object@list_of_col_data
		col_data <- data[names(data) %in% which_data]
		data_list <- list()
		for(i in 1:length(col_data)) {
			names <- col_data[[i]][ ,metab_ids[i]]
			groups <- col_data[[i]][ ,metab_groups[i]]
			class <- rep(which_data[i], each=length(names))
			data_list[[i]] <- as.data.frame(cbind(names, groups, class))
			colnames(data_list[[i]]) <- c("names", "groups", "class")
		}
		col_data <- do.call(rbind, data_list)
		col_data <- as.data.frame(col_data)
		colnames(new_data_for_pca) <- unlist(lapply(strsplit(colnames(new_data_for_pca), split="_data.", fixed=TRUE), function(x) return(x[2])))
		new_data_for_pca <- new_data_for_pca[ ,match(colnames(new_data_for_pca), col_data$names)]
		phenotype_name <- object@annotations$phenotype
		phenotype <- object@list_of_data[[phenotype_name]]
		phenotype <- phenotype[rownames(phenotype) %in% rownames(new_data_for_pca),]
		phenotype <- phenotype[match(rownames(phenotype), rownames(new_data_for_pca)),]
		data <- new_data_for_pca
		pca_metabs <- prcomp(t(data), scale.=T, center=T)
		pca_individuals <- prcomp(data, scale.=T, center=T)
		pca_metabs_pcs <- as.data.frame(pca_metabs$x[,1:2])
		plot_data_metabs <- as.data.frame(cbind(pca_metabs_pcs, col_data))
		colnames(plot_data_metabs) <- c("PC1", "PC2", "names","groups", "class")
		pca_plot_metabs <- ggplot(plot_data_metabs, aes(x=PC1, y=PC2, color=groups, shape=class)) + geom_point() + labs(x="PC1", y="PC2", subtitle="PCA plot for metabolites")
		pca_individuals_pc <- as.data.frame(pca_individuals$x[,1:2])
		timepoints <- unlist(lapply(strsplit(rownames(new_data_for_pca), split="_"), function(x) return(x[2])))
		levels <- unique(sort(as.numeric(unlist(lapply(strsplit(timepoints, split="t", fixed=TRUE), function(x) return(x[2]))))))
		plot_data_individuals <- as.data.frame(cbind(pca_individuals_pc, phenotype[, cols_for_vis]))
		plot_data_individuals$timepoints <- factor(timepoints, levels=paste("t", levels, sep=""))
		colnames(plot_data_individuals)[1:2] <- c("PC1","PC2")
		pca_plot_individuals <- plot_ly(plot_data_individuals, x=~PC1, y=~PC2, color=~timepoints, colors=get_palette(30)[c(16,2,3,6,7,8,10,13,19,20,23,27)], type = 'scatter', 
															mode = 'markers', hoverinfo = 'text', 
														text = get_text(data=plot_data_individuals, colnames=cols_for_vis))
		return(list(metabs=pca_plot_metabs, samples=pca_plot_individuals))
	})

#lol <- pca_plotter_wrt_common(object=data, which_data=c("lipid_data","nmr_data"), metab_groups=c("sub_pathway","Group"), metab_ids=c("metabolite","id"), cols_for_vis=cols_for_pca_plot)


#' Function to plot tSNE plots for multiple datasets
#' @param object An object of class metab_analyser
#' @param cols_for_vis choose the colnames from phenotype data that you want to show in the text labels
#' @param which_data a character vector - Names of the dataset from which the samples will be extracted
#' @param metab_ids a character vector - names of the column with metabolite names in the col data matix(please make sure they are in the same order as the above ones)
#' @param metab_groups a character vector - names of the column which has the pathway information of the metabolites(please make sure they are in the same order as the above ones)
#' @param phenotype_index index of the phenotype data. Can input either the name of the phenotype dataset or the index of the same
#' @return a list with two plot objects 1) samples - tSNE plot of the individuals(".$samples")
#'										 2) metabs - tSNE plot of the metabolites(".$metabs")
#' @export

setGeneric("tsne_plotter_wrt_common", function(object, which_data, metab_groups, metab_ids, cols_for_vis) standardGeneric("tsne_plotter_wrt_common"))

setMethod("tsne_plotter_wrt_common", "metab_analyser", function(object, which_data, metab_groups, metab_ids, cols_for_vis) {
	new_list_of_data <- common_sample_extractor(object)
	new_list_of_data <- new_list_of_data[names(new_list_of_data) %in% which_data]
	new_data_for_tsne <- do.call(cbind, new_list_of_data)
	new_data_for_tsne <- na.omit(new_data_for_tsne)
	new_data_for_tsne <- as.data.frame(new_data_for_tsne)
	data <- object@list_of_col_data
	col_data <- data[names(data) %in% which_data]
	data_list <- list()
	for(i in 1:length(col_data)) {
		names <- col_data[[i]][ ,metab_ids[i]]
		groups <- col_data[[i]][ ,metab_groups[i]]
		class <- rep(which_data[i], each=length(names))
		data_list[[i]] <- as.data.frame(cbind(names, groups, class))
		colnames(data_list[[i]]) <- c("names", "groups", "class")
	}
	col_data <- do.call(rbind, data_list)
	col_data <- as.data.frame(col_data)
	colnames(new_data_for_tsne) <- unlist(lapply(strsplit(colnames(new_data_for_tsne), split="_data.", fixed=TRUE), function(x) return(x[2])))
	new_data_for_tsne <- new_data_for_tsne[ ,match(colnames(new_data_for_tsne), col_data$names)]
	phenotype_name <- object@annotations$phenotype
	phenotype <- object@list_of_data[[phenotype_name]]
	phenotype <- phenotype[rownames(phenotype) %in% rownames(new_data_for_tsne),]
	phenotype <- phenotype[match(rownames(phenotype), rownames(new_data_for_tsne)),]
	data <- new_data_for_tsne
	tsne_fit_metabs <- tsne(data)
	plot_data_metabs <- as.data.frame(cbind(as.data.frame(tsne_fit_metabs$data), col_data))
	colnames(plot_data_metabs) <- c("X1", "X2", "names","groups", "class")
	tsne_plot_metabs <- ggplot(plot_data_metabs, aes(x=X1, y=X2, color=groups, shape=class)) + geom_point() + labs(x="X1", y="X2", subtitle="tSNE plot for metabolites")
	tsne_samples <- tsne(t(data))	
	timepoints <- unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[2])))
	levels <- unique(sort(as.numeric(unlist(lapply(strsplit(timepoints, split="t", fixed=TRUE), function(x) return(x[2]))))))
	tsne_data <- as.data.frame(cbind(as.data.frame(tsne_samples$data), phenotype[,cols_for_vis]))
	tsne_data$timepoints <- factor(timepoints, levels=paste("t", levels, sep=""))
	tsne_plot_samples <- plot_ly(tsne_data, x=~X1, y=~X2, color=~timepoints, colors=get_palette(30)[c(16,2,3,6,7,8,10,13,19,20,23,27)], type = 'scatter', 
															mode = 'markers', hoverinfo = 'text', 
														text = get_text(data=tsne_data, colnames=cols_for_vis))
	return(list(metabs=tsne_plot_metabs, samples=tsne_plot_samples))
})

#lol <- tsne_plotter_wrt_common(object=data, which_data=c("lipid_data","nmr_data"), metab_groups=c("sub_pathway","Group"), metab_ids=c("metabolite","id"), cols_for_vis=cols_for_pca_plot)

#' Function to plot UMAP plots for multiple datasets
#' @param object An object of class metab_analyser
#' @param cols_for_vis choose the colnames from phenotype data that you want to show in the text labels
#' @param which_data a character vector - Names of the dataset from which the samples will be extracted
#' @param metab_ids a character vector - names of the column with metabolite names in the col data matix(please make sure they are in the same order as the above ones)
#' @param metab_groups a character vector - names of the column which has the pathway information of the metabolites(please make sure they are in the same order as the above ones)
#' @param phenotype_index index of the phenotype data. Can input either the name of the phenotype dataset or the index of the same
#' @return a list with two plot objects 1) samples - UMAP plot of the individuals(".$samples")
#'										 2) metabs - UMAP plot of the metabolites(".$metabs")
#' @export

setGeneric("umap_plotter_wrt_common", function(object, which_data, metab_groups, metab_ids, cols_for_vis, phenotype_index) standardGeneric("umap_plotter_wrt_common"))

setMethod("umap_plotter_wrt_common", "metab_analyser", function(object, which_data, metab_groups, metab_ids, cols_for_vis, phenotype_index) {
	new_list_of_data <- common_sample_extractor(object)
	umap_data <- new_list_of_data[names(new_list_of_data) %in% which_data]
	umap_data <- do.call(cbind, umap_data)
	umap_data <- as.data.frame(na.omit(umap_data))
	colnames(umap_data) <- unlist(lapply(strsplit(colnames(umap_data), split="_data.", fixed=TRUE), function(x) return(x[2])))
	data <- object@list_of_col_data
	col_data <- data[names(data) %in% which_data]
	data_list <- list()
	for(i in 1:length(col_data)) {
		names <- col_data[[i]][ ,metab_ids[i]]
		groups <- col_data[[i]][ ,metab_groups[i]]
		class <- rep(which_data[i], each=length(names))
		data_list[[i]] <- as.data.frame(cbind(names, groups, class))
		colnames(data_list[[i]]) <- c("names", "groups", "class")
	}
	col_data <- do.call(rbind, data_list)
	col_data <- as.data.frame(col_data)
	umap_data <- umap_data[ ,match(colnames(umap_data), col_data$names)]
	phenotype_name <- object@annotations$phenotype
	phenotype <- object@list_of_data[[phenotype_name]]
	phenotype <- phenotype[rownames(phenotype) %in% rownames(umap_data),]
	phenotype <- phenotype[match(rownames(phenotype), rownames(umap_data)),]
	umap_fit_metabs <- umap::umap(t(umap_data))
	plot_data_metabs <- cbind(as.data.frame(umap_fit_metabs$layout), col_data)
	colnames(plot_data_metabs) <- c("UMAP1", "UMAP2", "names","groups", "class")
	umap_plot_metabs <- ggplot(plot_data_metabs, aes(x=UMAP1, y=UMAP2, color=groups, shape=class)) + geom_point() + labs(x="UMAP1", y="UMAP2", subtitle="UMAP plot for metabolites")
	umap_fit_samples <- umap::umap(umap_data)
	timepoints <- unlist(lapply(strsplit(rownames(umap_data), split="_"), function(x) return(x[2])))
	levels <- unique(sort(as.numeric(unlist(lapply(strsplit(timepoints, split="t", fixed=TRUE), function(x) return(x[2]))))))
	plot_data_samples <- as.data.frame(cbind(as.data.frame(umap_fit_samples$layout), phenotype[,cols_for_vis]))
	plot_data_samples$timepoints <- factor(timepoints, levels=paste("t", levels, sep=""))
	colnames(plot_data_samples)[1:2] <- c("UMAP1","UMAP2")
	umap_plot_samples <- plot_ly(plot_data_samples, x=~UMAP1, y=~UMAP2, color=~timepoints, colors=get_palette(30)[c(16,2,3,6,7,8,10,13,19,20,23,27)], type = 'scatter', 
															mode = 'markers', hoverinfo = 'text', 
														text = get_text(data=plot_data_samples, colnames=cols_for_vis))
	return(list(samples=umap_plot_samples, metabs=umap_plot_metabs))
})

#lol <- umap_plotter_wrt_common(object=object, which_data=c("lipid_data","nmr_data"), metab_groups=c("sub_pathway","Group"), metab_ids=c("metabolite","id"), cols_for_vis=cols_for_pca_plot)
