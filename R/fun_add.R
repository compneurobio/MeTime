
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


#' creating metime_plotter class that converts calculations and metadata as a plotable object to parse 
#' into viz_plotter
#' Contains slots - plot_data: Dataframe with plotting data and metadata for visualization
#' 				  - plot: ggplot(), circos() or visNetwork() object with predefined aesthetics 
#'                - calc_type: A vector to specify type of calculation - will be used for comp_ functions
#'                - calc_info: string to define the information about calculation
#' 				  - plot_type: A character vector to define the type of plots that are needed.
#' @rdname metime_plotter
#' @export
setClass("metime_plotter", slots=list(plot_data="list", plot="list", calc_type="character", calc_info="character", plot_type="character", style="character"))


#' Function to add measurements taken at screening time for samples to be added to all timepoints
#' @description A method applied on the s4 object of class "metime_analyser" to add all those datapoints that were measured only during screening
#' to all the respective samples at all timepoints
#' @examples # adding APOEGrp, PTGENDER to all data points
#' new_with_apoegrp_sex <- add_screening_vars(object=metime_analyser_object, vars=c("APOEGrp","PTGENDER"))
#' @param object An object of class metime_analyser
#' @param vars A character naming the vars of interest
#' @return phenotype data which can be replaced into the original object or use it separately with a different object
#' @export
setGeneric("add_screening_vars", function(object, vars) standardGeneric("add_screening_vars"))

setMethod("add_screening_vars", "metime_analyser", function(object, vars) {
	phenotype_name <- object@annotations$phenotype
	phenotype <- object@list_of_data[[phenotype_name]]
	screening <- phenotype[grep("-1|-2", rownames(phenotype)), ]
	new_rows <- as.data.frame(screening[, vars])
	new_rows <- na.omit(new_rows)
	new_rows <- new_rows[order(rownames(new_rows)), ]
	sample_names <- unlist(lapply(strsplit(rownames(new_rows), split="_"), function(x) return(x[1])))
	for(i in 1:length(sample_names)) {
		samples <- unlist(lapply(strsplit(rownames(phenotype), split="_"), function(x) return(x[1])))
		index <- which(samples %in% sample_names[i], arr.ind=TRUE)
		for(j in index) {
			phenotype[j, vars] <- new_rows[rownames(new_rows)[i], vars]
		}
	}
	object@list_of_data[[phenotype_name]] <- phenotype
	out <- object
	return(out)
})

#' Function to check normality and add data to col data
#' @description A method applied on the s4 object of class "metime_analyser" to check normality of the metabolites
#' and add it to corresponding columns
#' @examples 
#'	object <- add_col_normality(object=data, which_data=c("lipid_data","nmr_data"), type="shapiro", metab_names=c("metabolite","id"))
#' @param object An object of class metime_analyser
#' @param which_data dataset on which the method is to be applied
#' @param type type of test, "shapiro" and "kruskal" are available
#' @param metab_names column that has the metabolite names in col_data.
#' @param all logical to add all kinds of available stats.
#' @return S4 object with shapiro wilk test related data in the col_data
#' @export
setGeneric("add_col_stats", function(object, which_data, type, metab_names, all) standardGeneric("add_col_stats"))
setMethod("add_col_stats", "metime_analyser", function(object, which_data, type="shapiro", metab_names, all) {
	  stopifnot(type %in% c("shapiro", "kruskal"))
      for(i in 1:length(which_data)) {
      		out <- lapply(names(object@list_of_data[[which_data[i]]]), function(x) {
      					my_model <- shapiro.test(object@list_of_data[[which_data[i]]][, x])
      					return(data.frame(id=x,
                 				shapiro_pval=as.numeric(my_model$p.value),
                 				shapiro_statistic=as.numeric(my_model$statistic),
                 				stringsAsFactors = F))}) %>%  
      					do.call(what=rbind.data.frame)
	      	dummy <- out[order(out$id), ]
	      	shapiro_pval <- dummy$shapiro_pval
	      	shapiro_statistic <- dummy$shapiro_statistic
	      	shapiro_normal <- ifelse(dummy$shapiro_pval>0.05, TRUE,FALSE)
	      	out <- lapply(names(object@list_of_data[[which_data[i]]]), function(x) {
      						my_model <- shapiro.test(object@list_of_data[[which_data[i]]][ ,x])
      						return(data.frame(id=x,
                 				kruskal_pval=as.numeric(my_model$p.value),
                 				kruskal_statistic=as.numeric(my_model$statistic),
                 				stringsAsFactors = F))}) %>%  
      						do.call(what=rbind.data.frame)
      		dummy <- out[order(out$id), ]
	      	kruskal_pval <- dummy$kruskal_pval
	      	kruskal_statistic <- dummy$kruskal_statistic
	      	kruskal_normal <- ifelse(dummy$kruskal_pval>0.05, TRUE,FALSE)
	      	col_data <- object@list_of_col_data[[which_data[i]]]
	    	col_data <- col_data[order(col_data[ ,metab_names[i]]), ]
	    	out <- lapply(names(object@list_of_data[[which_data[i]]]), function(x) {
	    				data <- object@list_of_data[[which_data[i]]][ ,x]
	    				names(data) <- rownames(object@list_of_data[[which_data[i]]])
	    				times <- unlist(lapply(strsplit(names(data), split="_"), function(x) return(x[2])))
	    				out <- ifelse(length(unique(times))>1, TRUE, FALSE)
	    				return(out)
	    		})
	    	longitudinal <- unlist(out)
	    	object@list_of_col_data[[which_data[i]]]$longitudinal <- longitudinal
	      	if(all) {
	      		col_data <- as.data.frame(cbind(col_data, kruskal_pval, kruskal_statistic, kruskal_normal, 
	      										shapiro_pval, shapiro_statistic, shapiro_normal))
	      		object@list_of_col_data[[which_data[i]]] <- col_data
	      	} else {
	    		if(type %in% "shapiro") {
	    			col_data <- as.data.frame(cbind(col_data, shapiro_pval, shapiro_statistic, shapiro_normal))
	    			object@list_of_col_data[[which_data[i]]] <- col_data
	    		} else if(type %in% "kruskal") {
	    			col_data <- as.data.frame(cbind(col_data, kruskal_pval, kruskal_statistic, kruskal_normal))
	      			object@list_of_col_data[[which_data[i]]] <- col_data
	    		}  
	      	}
      		
    }  
    out <- object    	
	return(out)
  
})


#' Function to add measurements taken at screening time for samples to be added to all timepoints in row data
#' @description A method applied on the s4 object of class "metime_analyser" to add all those datapoints that were measured only during screening
#' to all the respective samples at all timepoints in row_data lists
#' @examples # adding APOEGrp, PTGENDER, and diag group to all data points and prepping the object for viz_distribution_plotter()
#' object <- add_distribution_vars_to_rows(object=data, screening_vars=c("APOEGrp", "DXGrp_longi", "PTGENDER"), 
#'			distribution_vars=c("Age", "BMI", "ADNI_MEM", "ADNI_LAN", "ADNI_EF", "APOEGrp", "DXGrp_longi", "PTGENDER"), which_data="lipid_data")
#' @param object An object of class metime_analyser
#' @param vars A character naming the vars of interest
#' @param which_data dataset to which the information is to be added(only 1 can be used at a time)
#' @return object of class metime_analyser with phenotype data added to row data
#' @export
setGeneric("add_distribution_vars_to_rows", function(object, screening_vars, distribution_vars, which_data) standardGeneric("add_distribution_vars_to_rows"))
setMethod("add_distribution_vars_to_rows", "metime_analyser", function(object, screening_vars, distribution_vars, which_data) {
			if(!is.null(screening_vars)) {
				phenotype <- add_screening_vars(object, screening_vars)
			} else {
				phenotype <- object@list_of_data[[object@annotations$phenotype]]
			}
			data <- as.data.frame(object@list_of_row_data[[which_data]])
			phenotype <- phenotype[rownames(phenotype) %in% rownames(data), ]
			phen_data <- phenotype[order(rownames(phenotype)),]
			data <- data[order(rownames(data)), ]
			data <- cbind(data, phen_data[ ,distribution_vars])
			object@list_of_row_data[[which_data]] <- data
			out <- object
			return(out)
	})

#' Function to add covariates to the dataset of interest for GGMs
#' @description adds Covariates to data matrices in metime_analyser S4 object
#' @param object object of class metime_analyser
#' @param which_data Dataset to which the covariates is to be added
#' @param covariates character vector names of covariates. 
#' @param class.ind Logical to convert factor variables into class.ind style or not
#' @param phenotype Logical. If True will extract from phenotype dataset else uses row data
#' @return S4 object with covariates added to the dataset
#' @export
setGeneric("add_phenotypes_as_covariates", function(object, which_data, covariates, class.ind, phenotype) standardGeneric("add_phenotypes_as_covariates"))
setMethod("add_phenotypes_as_covariates", "metime_analyser", function(object, which_data, covariates, class.ind, phenotype) {
			if(phenotype) {
				phenotype <- object@list_of_data[[object@annotations$phenotype]]
				phenotype <- phenotype[order(rownames(phenotype)),]
				list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
				covariates_list <- lapply(list_of_data, function(x) {
							x <- x[order(rownames(x)), ]
							return(phenotype[rownames(phenotype) %in% rownames(x), covariates])
					})
				for(i in 1:length(which_data)) {
					 list_of_data <- unname(list_of_data)
					 covariates_list <- unname(covariates_list)
					 list_of_data[[i]] <- as.data.frame(cbind(list_of_data[[i]], covariates_list[[i]]))
					 if(class.ind) {
					 	for(j in 1:length(covariates)) {
					 		vec <- list_of_data[[i]][ ,covariates[j]]
					 		if(class(vec) %in% "character" | class(vec) %in% "factor") {
					 			new_data <- class.ind(vec)
					 			list_of_data[[i]] <- list_of_data[[i]][ ,!(names(list_of_data[[i]]) %in% covariates[i])]
					 			list_of_data[[i]] <- as.data.frame(cbind(list_of_data[[i]], new_data))
							}
					 	}
					 } else {
					 	for(j in 1:length(covariates)) {
					 			vec <- list_of_data[[i]][ ,covariates[j]]
					 			if(class(vec) %in% "character" | class(vec) %in% "factor") {
					 				new_vec <- as.numeric(vec)
					 				list_of_data[[i]] <- list_of_data[[i]][ ,!(names(list_of_data[[i]]) %in% covariates[i])]
					 				list_of_data[[i]] <- as.data.frame(cbind(list_of_data[[i]], new_vec))
					 			}
					 		}
					 }
					 object@list_of_data[[which_data[i]]] <- list_of_data[[i]]
				}
			} else {
				list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
				for(i in 1:length(which_data)) {
					dummy <- list_of_row_data[[which_data[i]]]
					dummy <- dummy[ ,covariates[[i]]]
					dummy <- dummy[order(rownames(dummy)), ]
					data <- list_of_data[[i]]
					data <- data[rownames(data) %in% rownames(dummy), ]
					data <- data[order(rownames(data)), ]
					list_of_data[[i]] <- as.data.frame(cbind(data, dummy))
					if(class.ind) {
					 	for(j in 1:length(covariates)) {
					 		vec <- list_of_data[[i]][ ,covariates[j]]
					 		if(class(vec) %in% "character" | class(vec) %in% "factor") {
					 			new_data <- class.ind(vec)
					 			list_of_data[[i]] <- list_of_data[[i]][ ,!(names(list_of_data[[i]]) %in% covariates[i])]
					 			list_of_data[[i]] <- as.data.frame(cbind(list_of_data[[i]], new_data))
							}
					 	}
					 } else {
					 	for(j in 1:length(covariates)) {
					 			vec <- list_of_data[[i]][ ,covariates[j]]
					 			if(class(vec) %in% "character" | class(vec) %in% "factor") {
					 				new_vec <- as.numeric(vec)
					 				list_of_data[[i]] <- list_of_data[[i]][ ,!(names(list_of_data[[i]]) %in% covariates[i])]
					 				list_of_data[[i]] <- as.data.frame(cbind(list_of_data[[i]], new_vec))
					 			}
					 		}
					 }
					 object@list_of_data[[which_data[i]]] <- list_of_data[[i]]
				}
			}
			out <- object
			return(out)
	})	

#' Function to add metabolites as covariates for network construction
#' @description Method applied on metime_analyser object to add other metabolite data to a certain dataset
#' @param object A S4 object of class metime_analyser
#' @param which_data Dataset to which the metab data is to be added(please note that this a single character)
#' @param which_metabs list of names of metabs and name of the list represents the dataset from which 
#' the metabs are to be acquired. eg: which_metabs=list(nmr_data=c("metab1", "metab2"), lipid_data=c(""))
#' @return S4 object with metabs added for GGM to another dataset
#' @export
setGeneric("add_metabs_as_covariates", function(object, which_data, which_metabs) standardGeneric("add_metabs_as_covariates"))
setMethod("add_metabs_as_covariates", "metime_analyser", function(object, which_data, which_metabs) {
			data <- object@list_of_data[[which_data]]
			list_of_metabs <- list()
			for(i in 1:length(which_metabs)) {
				dat <- object@list_of_data[[names(which_metabs)[i]]]
				dat <- dat[ ,which_metabs[[i]]]
				dat <- dat[order(rownames(dat)), ]
				list_of_metabs[[i]] <- dat
			}
			if(length(which_data) > 1) {
				metab_matrix <- do.call(cbind, list_of_metabs)
			} else {
				metab_matrix <- list_of_metabs[[1]]
			}
			list_of_samples <- list(rownames(data), rownames(metab_matrix))
			common_samples <- Reduce(intersect, list_of_samples)
			data <- data[rownames(data) %in% common_samples, ]
			data <- data[order(rownames(data)),]
			metab_matrix <- metab_matrix[rownames(metab_matrix) %in% common_samples, ]
			metab_matrix <- metab_matrix[order(rownames(metab_matrix)), ]
			object@list_of_data[[which_data]] <- as.data.frame(cbind(data, metab_matrix))
			out <- object
			return(out)
	})


#' Function to add features to visnetwork plot from another plotter object
#' @description Function to add node features to see the nodes in the network that affected differently
#' @param network_plotter_object plotter object with network information
#' @param guide_plotter_object guide from which the colors are to be extracted
#' @param which_type type of the guide plotter object to be used. Current options are "regression" and "conservation"
#' @param metab_colname name of the column in guide plotter object that represents the metabolites
#' @return network plotter object with new node colors/features
#' @export
setGeneric("add_node_features", function(network_plotter_object, guide_plotter_object, which_type, metab_colname) standardGeneric("add_node_features"))
setMethod("add_node_features", "metime_plotter", function(network_plotter_object, guide_plotter_object, which_type, metab_colname) {
		network_plotter_object@plot_data[["node"]] <- network_plotter_object@plot_data[["node"]][order(network_plotter_object@plot_data[["node"]]$label), ]
		guide_plotter_object@plot_data[[1]] <- guide_plotter_object@plot_data[[1]][guide_plotter_object@plot_data[[1]][ ,metab_colname] %in% network_plotter_object@plot_data[[1]]$label, ]
		guide_plotter_object@plot_data[[1]] <- guide_plotter_object@plot_data[[1]][order(guide_plotter_object@plot_data[[1]][ ,metab_colname]), ]
		if(which_type %in% "regression") {
			column_for_colors <- guide_plotter_object@plot_data[[1]][ ,c("beta", "pval")]
			column_for_colors <- sign(column_for_colors$beta) * -log10(column_for_colors$pval)
		} else if(which_type %in% "conservation") {
			column_for_colors <- guide_plotter_object@plot_data[[1]][ ,"ci"]
		}
		color.gradient <- function(x, colors=c("blue","gray","red"), colsteps=50) {
  			return(colorRampPalette(colors) (colsteps)[findInterval(x, seq(-1, 1, length.out=colsteps))] )
		}
		gradient <- color.gradient(column_for_colors)
		network_plotter_object@plot_data[["node"]]$color <- gradient
		return(network_plotter_object)
	}) 