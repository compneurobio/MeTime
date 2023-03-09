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
      add <- lapply(seq_along(which_data), function(i) {
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
	      	return(NULL)
      		
    })  
    out <- object    	
	return(out)
  
})


