usethis::use_import_from("rstatix", fun=c("wilcox_test", "t_test"))
usethis::use_package("tidyverse", type="depends")
usethis::use_package("DT", type="Imports")
usethis::use_import_from("M3C", fun="tsne") # Also has umap - will change it
usethis::use_import_from("umap", fun="umap")
usethis::use_import_from("plotly", fun="ggplotly")
usethis::use_package("visNetwork", type="Imports")
usethis::use_package("GeneNet", type="Imports")
usethis::use_package("longitudinal", type="Imports")
usethis::use_package("glmnet", type="Imports")
usethis::use_package("Boruta", type="Imports")
usethis::use_import_from("multiway", fun="parafac") # get that manually - can't dependence
usethis::use_import_from("xlsx", fun="write.xlsx") # No base function and same dependence issue
usethis::use_package("mgcv", type="Imports")
usethis::use_import_from("lmerTest", fun="lmer")
usethis::use_import_from("parallel", fun="mclapply")
usethis::use_package("WGCNA", type="Imports") # can't dependence
usethis::use_import_from("dynamicTreeCut", fun="cutreeDynamic")
usethis::use_package("rmarkdown", type="Imports")
usethis::use_package("rmdformats", type="Imports")
usethis::use_package("htmltools", type="Imports")
usethis::use_package("rlang", type="Imports")
usethis::use_package("MatrixEQTL", type="Imports")


#' Validity function to check if the object is valid or not
#' @description Function to check the validity of metime_analyser 
#' @param object An S4 object of class metime_analyser
#' @examples validObject(object)
#' @returns logical suggesting if the object is intact or not
#' @export
validity <- function(object) {
		out <- TRUE
		if(!check_rownames_and_colnames(object)) out <- FALSE 
		if(!check_ids_and_classes(object)) out <- FALSE
		if(!check_results(object)) out <- FALSE
		return(out) 
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
								 annotations="list", results="list"), validity=validity) 



# add a function to generate Rmarkdown directly from the analyser object - viz_analyser()
# create another function which add the function applied on the analyser_object - ongoing
# Try to improve the code of all functions and then Matthias will check it

#' Setting a plotting method for the metime_analyser class
#' @description Function to plot results of a certain calculation 
#' @param x An S4 object of class metime_analyser
#' @param results_index Index/name of the results to be plotted
#' @param interactive logical. Set TRUE for interactive plot
#' @param ... other parameters to pass color, fill, strat, viz(character vector with colnames for interactive)
#' @param plot_type to define the type of plot. Accepted inputs are "dot", "tile", "box", "forrest", "manhattan"
#' @return plots for a certain set of results
#' @export
setMethod("plot", "metime_analyser", function(x, results_index, interactive, plot_type, ...) {
		add <- list(...)
		results <- x@results[[results_index]]
		get_text_for_plot <- function(data, colnames) {	
						out <- c()
						count <- 1
						text <- c()
						for(l in 1:length(rownames(data))) {
							for(m in 1:length(colnames)) {
								if(m==1) {
									text <- paste("<br /> ", colnames[m], " : ", data[l, colnames[m]], sep="")
								} else {
									text <- paste(text, "<br /> ", colnames[m], " : ", data[l, colnames[m]], sep="")
								}
							}
							out[count] <- text
							count <- count + 1
						} 
						return(out)
			} 

		make_interactive <- function(.plot, j, col, .results=results) {
					plot <- plotly::ggplotly(.plot, width=800, height=800)
					for(i in seq_along(plot$x$data)) {
						if(!is.null(plot$x$data[[i]]$text)) {
							x <- plot$x$data[[i]]$x
							data <- .results$plot_data[[j]]
							data <- data[data[ ,col] %in% x, ] 
							plot$x$data[[i]]$text <- get_text_for_plot(data=data, 
										colnames=colnames(data))
						}
					}
					return(plot)
			}

		if(is.null(names(results$plot_data))) {
			all_plots <- lapply(seq_along(results$plot_data), function(ind_data) {
					if(results$information$calc_type[ind_data] %in% "CI_metabolite" |
						results$information$calc_type[ind_data] %in% "CI_metabotype") {
						df <- results$plot_data[[ind_data]]
						exclude_cols <- c("x","y", "ci", "id", "time_from", "time_to", "nsubject", "n", "rank", "cor", "id_from", "id_to", "samples", "timepoints")
						df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.character) | sapply(df, is.factor))] <- apply(df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.character) | sapply(df, is.factor))], 2, factor)
						df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.numeric) | sapply(df, is.integer))] <- lapply(df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.numeric) | sapply(df, is.integer))], function(a) {
  										if(is.na(min(a)) | is.na(max(a))) {
  											breaks <- c(-Inf, quantile(a, probs = c(0.25, 0.5, 0.75), na.rm=TRUE), Inf)
  										} else {
  											breaks <- c(min(a), quantile(a, probs = c(0.25, 0.5, 0.75), na.rm=TRUE), max(a))
  										}
  										return(cut(a, breaks = breaks, labels = c("Q1", "Q2", "Q3", "Q4")))
						})
						if(plot_type %in% "dot") {
							plot <- ggplot(df, aes(x=x, y=ci)) + 
								geom_point(aes_string(color=add$color, 
								shape=add$shape)) + facet_wrap(add$strats) + theme_classic() +
								scale_shape_manual(values=1:nlevels(as.factor(df[ ,add$shape]))) + xlab("Subject order")
						} else if(plot_type %in% "box") {
							if(is.null(add$box_x)) {
								stop("Define x-axis for boxplot. Use box_x=variable as an argument in the plot function")
							} 
							median_data <- df
							median_data <- median_data[median_data$ci != 1, ]
							median_data$ci <- -log10(1-median_data$ci)
							combinations <- combn(unique(median_data[ ,add$box_x]), 2)
							if(length(which(is.na(df[ ,add$box_x])))>=1) {
								info <- paste("There are ", length(which(is.na(df[ ,add$box_x]))), "NAs in this variable and are removed")
							} else {
								info <- NULL
							}
							median_data <- median_data[!is.na(median_data[ ,add$box_x]), ]	
							my_comparisons <- lapply(1:ncol(combinations), function(x) return(as.character(combinations[,x])))
							plot <- ggpubr::ggboxplot(median_data, x=add$box_x, 
								y="ci", color="black", fill=add$box_x, alpha=0.5, notch=F) +
          						ggpubr::stat_compare_means(method="wilcox.test", comparisons=my_comparisons, 
          						label.y=seq(from=3.5, to=11, by=0.5)) 
          					plot <-	ggpubr::ggpar(plot, ylab="-log10(1-ci)") + labs(caption=info)
						} else {
							stop("This type of plot is not available for this calculation")
						}
						if(interactive) {
							plot <- make_interactive(.plot=plot, j=ind_data, col="x")
						} 
						return(plot)
					} else if(results$information$calc_type[ind_data] %in% "pairwise_distance" |
						results$information$calc_type[ind_data] %in% "pairwise_correlation") {
						if(plot_type %in% "tile") {
							plot <- ggplot(results$plot_data[[ind_data]], aes(x=row,y=column)) + 
							geom_tile(aes(fill=dist)) +
								facet_wrap(add$strats) + theme_classic()
						} else {
							stop("This type of plot is not available for this calculation")
						}
						if(interactive) {
							plot <- make_interactive(.plot=plot, j=ind_data, col="row")
						} 
						return(plot)
					} else if(results$information$calc_type[ind_data] %in% "colinearity") {
						if(plot_type %in% "tile") {
							plot <- ggplot(results$plot_data[[ind_data]], aes(x=med_class_1,y=med_class_2)) + 
							geom_tile(aes(fill=Cramers_V)) +
								facet_wrap(add$strats) + theme_classic()
						} else {
							stop("This type of plot is not available for this calculation")
						}
					} else if(results$information$calc_type[ind_data] %in% "PCA") {
						df <- results$plot_data[[ind_data]]
						exclude_cols <- c("PC1", "PC2", "time", "id", "subject", "name")
						df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.character) | sapply(df, is.factor))] <- apply(df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.character) | sapply(df, is.factor))], 2, factor)
						df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.numeric) | sapply(df, is.integer))] <- lapply(df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.numeric) | sapply(df, is.integer))], function(a) {
  										if(is.na(min(a)) | is.na(max(a))) {
  											breaks <- c(-Inf, quantile(a, probs = c(0.25, 0.5, 0.75), na.rm=TRUE), Inf)
  										} else {
  											breaks <- c(min(a), quantile(a, probs = c(0.25, 0.5, 0.75), na.rm=TRUE), max(a))
  										}
  										return(cut(a, breaks = breaks, labels = c("Q1", "Q2", "Q3", "Q4")))
						})
						if(plot_type %in% "dot") {
							plot <- ggplot(df, aes(x=PC1, y=PC2)) + 
								geom_point(aes_string(color=add$color, 
								shape=add$shape)) + facet_wrap(add$strats) + theme_classic() +
								scale_shape_manual(values=1:nlevels(as.factor(df[ ,add$shape])))
						} else {
							stop("This type of plot is not available for this calculation")
						}
						if(interactive) {
							plot <- make_interactive(.plot=plot, j=ind_data, col="PC1")
						} 
						return(plot)
					} else if(results$information$calc_type[ind_data] %in% "UMAP") {
						df <- results$plot_data[[ind_data]]
						exclude_cols <- c("UMAP1", "UMAP2", "time", "id", "subject", "name")
						df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.character) | sapply(df, is.factor))] <- apply(df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.character) | sapply(df, is.factor))], 2, factor)
						df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.numeric) | sapply(df, is.integer))] <- lapply(df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.numeric) | sapply(df, is.integer))], function(a) {
  										if(is.na(min(a)) | is.na(max(a))) {
  											breaks <- c(-Inf, quantile(a, probs = c(0.25, 0.5, 0.75), na.rm=TRUE), Inf)
  										} else {
  											breaks <- c(min(a), quantile(a, probs = c(0.25, 0.5, 0.75), na.rm=TRUE), max(a))
  										}
  										return(cut(a, breaks = breaks, labels = c("Q1", "Q2", "Q3", "Q4")))

						})
						if(plot_type %in% "dot") {
							plot <- ggplot(df, aes(x=UMAP1, y=UMAP2)) + 
								geom_point(aes_string(color=add$color, 
								shape=add$shape)) + facet_wrap(add$strats) + theme_classic() +
								scale_shape_manual(values=1:nlevels(as.factor(df[ ,add$shape])))
						} else {
							stop("This type of plot is not available for this calculation")
						}
						if(interactive) {
							plot <- make_interactive(.plot=plot, j=ind_data, col="UMAP1")
						} 
						return(plot)
					} else if(results$information$calc_type[ind_data] %in% "tSNE") {
						df <- results$plot_data[[ind_data]]
						exclude_cols <- c( "X1", "X2", "time", "id", "subject", "name")
						df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.character) | sapply(df, is.factor))] <- apply(df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.character) | sapply(df, is.factor))], 2, factor)
						df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.numeric) | sapply(df, is.integer))] <- lapply(df[, !(colnames(df) %in% exclude_cols) & (sapply(df, is.numeric) | sapply(df, is.integer))], function(a) {
  										if(is.na(min(a)) | is.na(max(a))) {
  											breaks <- c(-Inf, quantile(a, probs = c(0.25, 0.5, 0.75), na.rm=TRUE), Inf)
  										} else {
  											breaks <- c(min(a), quantile(a, probs = c(0.25, 0.5, 0.75), na.rm=TRUE), max(a))
  										}
  										return(cut(a, breaks = breaks, labels = c("Q1", "Q2", "Q3", "Q4")))
						})
						if(plot_type %in% "dot") {
							plot <- ggplot(df, aes(x=X1, y=X2)) + 
								geom_point(aes_string(color=add$color, 
								shape=add$shape)) + facet_wrap(add$strats) + theme_classic()+
								scale_shape_manual(values=1:nlevels(as.factor(df[ ,add$shape])))
						} else {
							stop("This type of plot is not available for this calculation")
						}
						if(interactive) {
							plot <- make_interactive(.plot=plot, j=ind_data, col="X1")
						} 
						return(plot)
					} else if(results$information$calc_type[ind_data] %in% "regression") {
						if(is.null(add$group)) add$group <- "y"
						if(plot_type %in% "forest") {
							if(!all(c("xmin", "xmax") %in% colnames(results$plot_data[[ind_data]]))) {
								plot <- ggplot(results$plot_data[[ind_data]], 
									aes_string(x="x", y=add$group, color="color", 
										fill="color")) +
        								geom_point() + scale_color_manual(name="color", 
        									values=c("none"="#EAE4E3","nominal"="#FCF6A4","li"="#D4F582","fdr"="#82DEF5","bonferroni"="#EE6868"))
							} else {
								plot <- ggplot(results$plot_data[[ind_data]], 
										aes_string(x="x", y=add$group, xmin="xmin", xmax="xmax", 
										color="color", fill="color")) +
        								geom_point() +
        								geom_errorbar(width=.2, position=position_dodge(0.05))  + scale_color_manual(name="color", 
        									values=c("none"="#EAE4E3","nominal"="#FCF6A4","li"="#D4F582","fdr"="#82DEF5","bonferroni"="#EE6868"))
							}
						} else {
							stop("This type of plot is not available for this calculation")
						}
						if(interactive) {
							plot <- make_interactive(.plot=plot, j=ind_data, col="x")
						}
						return(plot) 
					} else if(results$information$calc_type[ind_data] %in% "distribution_metabs") {
						data <- results$plot_data[[ind_data]]
						if(is.null(add$col)) {
							stop("col is not specified for distribution plot. Please set col=variable in the plot function")
						}
						if(plot_type %in% "bar") {
							plot_data <- data[ ,add$col]
							plot_data <- table(plot_data) %>% as.data.frame()
							colnames(plot_data) <- c(add$col, "Frequency")
							bar_plot <- ggplot(plot_data, aes_string(x=add$col, y="Frequency")) +
												geom_bar(stat="identity") + theme_classic() + coord_flip()
							return(bar_plot)
						} else if(plot_type %in% "density") {
							plot_data <- data[ ,add$col]
							mean <- mean(data[,add$col])
							density_plot <- ggplot(plot_data, aes_string(x=add$col)) + geom_density(alpha=0.5) +
											geom_vline(data=plot_data, aes(xintercept=mean), linetype="dashed") + theme_classic() +
											ylab("Density")
							return(density_plot) 
						} else {
							stop("This type of plot is not available for this calculation")
						}
					} else if(results$information$calc_type[ind_data] %in% "distribution_samples") {
						data <- results$plot_data[[ind_data]]
						if(is.null(add$col)) {
							stop("col is not specified for distribution plot. Please set col=variable in the plot function")
						}
						if(plot_type %in% "bar") {
							plot_data <- data[ ,c(add$col, "time")]
							plot_data <- table(plot_data)
							plot_data <- reshape2::melt(plot_data)
							colnames(plot_data) <- c(add$col, "Timepoints", "Frequency")
							levels <- plot_data$Timepoints %>% unique() %>% gsub(pattern="[a-z|A-Z]", replacement="") %>% as.numeric() %>% sort()
							plot_data$Timepoints <- factor(plot_data$Timepoints, levels=paste("t", levels, sep=""))
							bar_plot <- ggplot(data=plot_data, aes_string(x=add$col, y="Frequency", fill="Timepoints")) +
											geom_bar(stat="identity") + theme_classic()
							return(bar_plot)
						} else if(plot_type %in% "density") {
							plot_data <- data[ ,c(add$col, "time")]
							levels <- data$time %>% unique() %>% gsub(pattern="[a-z|A-Z]", replacement="") %>% as.numeric() %>% sort()
							plot_data$time <- factor(plot_data$time, levels=paste("t", levels, sep=""))
							colnames(plot_data) <- c(add$col, "Timepoints")
							plot_data <- na.omit(plot_data)
							mu <- as.data.frame(aggregate(plot_data[,add$col], list(Timepoints=plot_data$Timepoints), FUN=mean))
							colnames(mu)[2] <- "mean"
							density_plot <- ggplot(plot_data, aes_string(x=add$col, fill="Timepoints")) + geom_density(alpha=0.3) + 
										geom_vline(data=mu, aes(xintercept=mean, color=Timepoints), linetype="dashed") +
			 							theme_classic() + ylab("Density")
			 				return(density_plot)
						} else {
							stop("This type of plot is not available for this calculation")
						}
					} else if(results$information$calc_type[ind_data] %in% "feature_selection") {
						data <- results$plot_data[[ind_data]]
						if(plot_type %in% "manhattan") {
							
						}
					}
				})
				names(all_plots) <- results$information$calc_info
				return(all_plots)
			} else {
				node_list <- results$plot_data$network$node
				edge_list <- results$plot_data$network$edge
				#Choosing colors and shapes for visualization
				if(length(unique(node_list$class))>1) {
					shapes <- c("square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond")
					shapes_data <- data.frame(name=unique(node_list$class), shape=shapes[1:length(unique(node_list$class))])
					node_list <- node_list %>%
						dplyr::mutate(shape = shapes_data$shape[match(class, shapes_data$name)])				
				}

        		if(is.null(node_list$color)) {
        			if(!is.null(add$color)) {
        				column <- add$color
        			} else {
        				column <- "class"
        			}
            		if(is.numeric(node_list[[column]])){
    					# Create a gradient of colors using colorRampPalette
    					n_colors <- length(unique(node_list[[column]]))
    					colors <- colorRampPalette(rainbow(n_colors))(n_colors)
    					node_list$color <- colors
  					} else {
    					# Creating a contrasting color palette for pathway coloring
    					unique_vals <- unique(node_list[[column]])
    					num_vals <- length(unique_vals)
    					hues <- seq(15, 375, length.out = num_vals)
    					chroma <- 90
    					luminance <- 65
    					colors <- grDevices::hcl(h = hues, c = chroma, l = luminance)
    					names(colors) <- unique_vals
    					node_list$color <- colors[as.character(node_list[[column]])]	
  					}
        			
        		} else {
        			if(is.null(add$group)) {
        				column=NULL
        			} else {
        				column=add$group
        			}
        		}
        		graph <- visNetwork::visNetwork(nodes=node_list, edges=edge_list) %>%
                    	visNetwork::visIgraphLayout(layout="layout.fruchterman.reingold", physics = F, smooth = F) %>%
                    	visNetwork::visPhysics(stabilization = FALSE) %>%
                    	visNetwork::visLegend(useGroups = T) %>% 
                    	visNetwork::visNodes(borderWidth = 3) %>%
                    	visNetwork::visEdges(smooth = FALSE, shadow = TRUE) %>%
                    	visNetwork::visOptions(highlightNearest = list(enabled=T, hover=T), nodesIdSelection = T, selectedBy=column) %>%
                    	visNetwork::visInteraction(navigationButtons = T) %>%
                    	visNetwork::visExport(type="pdf", name=ifelse(is.null(add$title), "network_metime", add$title), 
                    		float="right")	
                    	return(list(network=graph))
        	}
        	
	})

#' Setting new print definition for the metime_analyser object
#' @description function to see the structure of metime_analyser object
#' @param object S4 object of class metime_analyser
#' @examples structure(object)
#' @return structure of the S4 object
#' @export
setMethod("show", "metime_analyser", function(object) {
		cat("Datasets: \n")
		out <- lapply(seq_along(object@list_of_data), function(i) {
			if(names(object@list_of_data)[i] %in% object@annotations$phenotype) {
				cat(paste(" - Dataset:", names(object@list_of_data)[i], "with", 
				dim(object@list_of_data[[i]])[1], "samples", "and", dim(object@list_of_data[[i]])[2], "phenotypes",sep=" "))
				cat("\n")
			} else if(names(object@list_of_data)[i] %in% object@annotations$medication) {
				cat(paste(" - Dataset:", names(object@list_of_data)[i], "with", 
				dim(object@list_of_data[[i]])[1], "samples", "and", dim(object@list_of_data[[i]])[2], "medicines", sep=" "))
				cat("\n")
			} else {
				cat(paste(" - Dataset:", names(object@list_of_data)[i], "with", 
				dim(object@list_of_data[[i]])[1], "samples", "and", dim(object@list_of_data[[i]])[2], "metabolites",sep=" "))
				cat("\n")
			}
		})
		cat("Results: \n")
		if(is.null(names(object@results))) return()
		out2 <- lapply(seq_along(object@results), function(j) {
				if(names(object@results)[j] %in% "") {
					return()
				}
				cat(paste(j, names(object@results)[j], sep="."))
				cat(" : ")
				cat("\n")
				cat(paste(" - ", object@results[[j]]$information$calc_info, collapse="\n"), sep="")
				cat("\n")
			})
	})


#' Function to stratify before calculation
#' @description Stratification of a dataset in calculations. Used only in calc_ functions
#' @param object An S4 object of class metime_analyser
#' @param which_data Dataset to be used
#' @param stratifications list of variables and their values to stratified - list(name=value)
#' @returns data with stratifications aforementioned
#' @export
setGeneric("get_stratified_data", function(object, which_data, stratifications) standardGeneric("get_stratified_data"))
setMethod("get_stratified_data", "metime_analyser", function(object, which_data, stratifications) {
		if(length(which_data)==1) {
      		data <- object@list_of_data[[which_data]]
      		row_data <- object@list_of_row_data[[which_data]]
      		col_data <- object@list_of_col_data[[which_data]]
		} else {
			object <- mod_extract_common_samples(object)
			object@list_of_data <- lapply(object@list_of_data, function(x) return(x[order(rownames(x)), ]))
			list_of_data <- object@list_of_data[which_data]
			list_of_col_data <- object@list_of_col_data[which_data]
			data <- do.call(cbind, unname(list_of_data))
			col_data <- do.call(dplyr::rbind.fill, unname(list_of_col_data))
			row_data <- object@list_of_row_data[[which_data[1]]]
			row_data <- row_data[rownames(row_data) %in% rownames(data), ]
		}
		if(length(stratifications)>=1) {
			row_data <- lapply(names(stratifications), function(x) {
              		row_data <- row_data[row_data[,x] %in% stratifications[[x]], ]
              		return(row_data) 
          		}) %>% do.call(what=rbind.data.frame)
        	data <- data[rownames(data) %in% rownames(row_data), ]		
		}
		return(list(data=data, row_data=row_data, col_data=col_data))
	})

#' Function to update plots post calculations
#' @description Function to update plots based on the calculation. 
#' See a calc_* function to see the usage of this function
#' @param object An S4 object of class metime_analyser
#' @param .interactive logical to make the plot interactive or not
#' @param results_index character/numeric input to define the results that you want to plot or replot with our automation.
#' Default will be set to NULL
#' @param type character to define the type of calculation used for updating the plot
#' cols_for_samples, cols_for_metabs, cols_for_meta etc will be used. So make sure you set those correctly
#' for better results. 
#' Allowed inputs are: c("ggm|network", "dimensionality_reduction", "CI_metabotype", "CI_metabolite", 
#' "pairwise_distance", "pairwise_correlation", "colinearity", "distribution", "regression", "feature_selection")
#' @returns object with plots of the newest calculation
#' @export

setGeneric("update_plots", function(object, .interactive=FALSE, type, results_index=NULL) standardGeneric("update_plots"))
setMethod("update_plots", "metime_analyser", function(object, .interactive=FALSE, type, results_index=NULL) {
		
		if(is.null(results_index)) {
			results <- object@results[[length(object@results)]]
		} else {
			results <- object@results[[results_index]]
		}

		if(grep("ggm|network", type) %>% length() == 1) {
			# chooses colors as combinations
			cols_of_int <- colnames(results$plot_data$network$node)
			cols_of_int <- cols_of_int[!cols_of_int %in% c("id", "label", "title")]
			if(!"color" %in% cols_of_int) {
				networks <- lapply(cols_of_int, function(a) {	
					network <- plot(object, results_index=length(object@results), 
						interactive=.interactive, plot_type="network", color=a)
					return(network)
				})
				names(networks) <- cols_of_int
			} else {
				networks <- list(node_feature=plot(object, results_index=length(object@results), 
						interactive=.interactive, plot_type="network"))
			}
			if(length(results$plots)==0) {
				results$plots[[1]] <- networks
			} else {
				results$plots[[length(results$plots)+1]] <- networks
			}
			object@results[[length(object@results)]] <- results
		} else if(type %in% "dimensionality_reduction") {
			cols_of_int <- colnames(results$plot_data[[1]])[!colnames(results$plot_data[[1]]) %in%
			c("PC1", "PC2", "UMAP1", "UMAP2", "X1", "X2", "time", "id", "subject", "name")]
			if(length(cols_of_int) < 2) {
				if(length(cols_of_int)==1) {
					combinations <- matrix(c(cols_of_int, cols_of_int), nrow=2, ncol=1)
				} else if(length(cols_of_int)==0) {
					plots <- list(no_meta=plot(object, results_index=length(object@results), 
						interactive=.interactive, plot_type="dot")) 
				}
			} else {
				combinations <- combn(cols_of_int, 2)
			}
			plots <- list()
			plots <- lapply(1:ncol(combinations), function(col) {
					dr_plots <- plot(object, results_index=length(object@results), 
						interactive=.interactive,
						plot_type="dot",
						color=as.character(combinations[1, col]),
						shape=as.character(combinations[2, col]))
					return(dr_plots)
				})
			names(plots) <- apply(combinations, 2, paste, collapse="-")
			if(length(results$plots)==0) {
				results$plots[[1]] <- plots
			} else {
				results$plots[[length(results$plots)+1]] <- plots
			}
			object@results[[length(object@results)]] <- results
		} else if(type %in% "CI_metabotype" | type %in% "CI_metabolite") {
			cols_of_int <- colnames(results$plot_data[[1]])[!colnames(results$plot_data[[1]]) %in% c("x","y", "ci", "id", "time_from", "time_to", "nsubject", "n", "rank", "cor", "id_from", "id_to", "samples", "timepoints", "name")]
			if(length(cols_of_int) < 2) {
				if(length(cols_of_int)==1) {
					combinations <- matrix(c(cols_of_int, cols_of_int), nrow=2, ncol=1)
				} else if(length(cols_of_int)==0) {
					plots <- list(no_meta=plot(object, results_index=length(object@results), 
						interactive=.interactive, plot_type="dot"))
					results$plots[[1]] <- plots 
				}
			} else {
				combinations <- combn(cols_of_int, 2)
			}
			dotplots <- lapply(1:ncol(combinations), function(col) {
					dot_plots <- plot(object, results_index=length(object@results), 
						interactive=.interactive,
						plot_type="dot",
						color=as.character(combinations[1, col]),
						shape=as.character(combinations[2, col]))
					return(dot_plots)
				})
			names(dotplots) <- apply(combinations, 2, paste, collapse="-")
			cols_of_int <- vapply(cols_of_int, function(col) {
						if(length(unique(results$plot_data[[1]][ ,col]))!=1) {
							return(col)
						} else {
							return("No")
						}
				}, character(1))
			cols_of_int <- cols_of_int[!cols_of_int %in% "No"]
			boxplots <- lapply(cols_of_int, function(col) {
						box_plots <- plot(object, results_index=length(object@results), 
							interactive=.interactive,
							plot_type="box",
							box_x=col)
						return(box_plots)
				})
			names(boxplots) <- cols_of_int
			if(length(results$plots)==0) {
				results$plots[[1]] <- dotplots
				results$plots[[2]] <- boxplots
			} else {
				results$plots[[length(results$plots)+1]] <- dotplots
				results$plots[[length(results$plots)+1]] <- boxplots
			}
			object@results[[length(object@results)]] <- results
		} else if(type %in% "pairwise_distance" | type %in% "pairwise_correlation" | type %in% "colinearity") {
			plots <- plot(object, results_index=length(object@results), interactive=.interactive,
				plot_type="tile")
			plots <- list(plots)
			names(plots) <- results$information$calc_info
			if(results$plots %>% length() == 1) {
				results$plots[[1]] <- plots
			} else {
				results$plots[[length(results$plots)+1]] <- plots
			}
			object@results[[length(object@results)]] <- results
		} else if(type %in% "regression") {
			data <- results$plot_data[[1]]
			cols_of_int <- colnames(data)[!colnames(data) %in% 
				c("x", "y", "xmin", "xmax", "pval", "tval", 
					"beta", "trait", "met", "se", "level", 
					"statistic","pvalue","FDR", "time", "type", "color")]
			if(length(cols_of_int) >= 1) {
				plots <- lapply(cols_of_int, function(x) {
						plots <-  plot(object, results_index=length(object@results), 
								interactive=.interactive, plot_type="forest", group=x)
						return(plots)
					})
				names(plots) <- cols_of_int
			} else {
				plots <-  plot(object, results_index=length(object@results), 
					interactive=.interactive, plot_type="forest")
				plots <- list(no_meta=plots)
			}
			if(results$plots %>% length() == 1) {
				results$plots[[1]] <- plots
			} else {
				results$plots[[length(results$plots)+1]] <- plots
			}
			object@results[[length(object@results)]] <- results
		} else if(type %in% "distribution") {
			data <- results$plot_data[[1]]
			plots <- lapply(colnames(data)[!colnames(data) %in% c("id", "time", "subject")], function(.col) {
						vec <- data[,.col] 
						vec <- na.omit(vec)
						if(is.character(vec)==TRUE | is.factor(vec)==TRUE) var_type <- "bar"
						if(is.numeric(vec)==TRUE) {
							if(length(unique(vec)) <= 20) {
								var_type <- "bar"
							} else {
								var_type <- "density"
							}
						}
						plots <- plot(object, interactive=.interactive, 
							results_index=length(object@results), col=.col, plot_type=var_type)
						return(plots)
					})
			names(plots) <- colnames(data)[!colnames(data) %in% c("id", "time", "subject")]
			if(results$plots %>% length() == 0) {
				results$plots[[1]] <- plots
			} else {
				results$plots[[length(results$plots)+1]] <- plots
			}
			object@results[[length(object@results)]] <- results
		} else if(type %in% "feature_selection") {
			data <- results$plot_data
		}
		return(object)
	})


#### ADD quantiles for continuous variables - done	

### GAMMs and LMMs to the code - done
### Plots for GAMMS - code ready to test
### add_data function into the package - done
### Manhattan plot into the package - need to do 
### testing the other networks - need to do
### 
