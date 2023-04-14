#' comparison of network analysis
#' @description A function to compare results from network analysis (GGMs)
#' @param object An S4 object of class metime_analyser
#' @param results_index character vector or numeric vector of length at least 2 which contains indices of results of interest
#' @param method character vector to define methods of interest. Accepted inputs are: "sign" for sign test and "het" for 
#' heterogenity test
#' @param name a character containing the name of the new result element appended to the analyzer object. 
#' Default set to 'meta_network_1'.
#' @return object with meta analyses of networks appended to results
#' @seealso [meta_parametric_test], [meta_nonparametric_test], [meta_regression]
#' @export
setGeneric("meta_network", function(object, results_index, method, name="meta_network_1") standardGeneric("meta_network"))
setMethod("meta_network", "metime_analyser", function(object, results_index, method, name="meta_network_1") {
		if(!all(method %in% c("sign", "cor","het"))) {
    		warning('method has to be one of "sign" (sign of beta) and "het" (heterogeneity). Exiting without making any changes')
    		return(object)
    	}
    	if(grep(name, names(object@results)) %>% length() == 1) {
    		warning("name of the results was previously used, using a different name")
      		index <- name %>% gsub(pattern="meta_network_", replacement="") %>% as.numeric() +1
      		name <- paste0("meta_network_",index)
    	}
    	if(is.null(results_index)){
      		results_index <- lapply(names(object@results), function(x) {
      				if(grep("ggm|network", object@results[[x]]$information$calc_type)==1) return(x)
      			}) %>% unlist()
    	}
    	my_data <- lapply(object@results[results_index], function(result) {
    			ret <- result$plot_data$network$edge
    			uids <- c()
    			for(i in 1:length(ret$pcor)) {
    				uids[i] <- paste(stringr::str_sort(c(ret$node1[i], ret$node2[i])), collapse="_")
    			}
    			ret$uids <- uids
    			return(ret) 
    		})
    	my_combn <- combn(names(my_data), 2) %>% t() %>% as.data.frame() %>% setNames(c("result1","result2"))
    	if("sign" %in% method) {
    		this_out <- lapply(1:nrow(my_combn), function(x){
        		this_result <- full_join(my_data[[my_combn$result1[x]]] %>%  # join data by met column
                                 dplyr::select(uids, pcor) %>% 
                                 dplyr::rename(uids=uids, pcor1=pcor) %>% 
                                 dplyr::mutate(sign1=ifelse(pcor1>=0, "+","-")),
                               my_data[[my_combn$result2[x]]] %>% 
                                 dplyr::select(uids, pcor) %>% 
                                 dplyr::rename(uids=uids, pcor2=pcor)%>% 
                                 dplyr::mutate(sign2=ifelse(pcor2>=0, "+","-")),
                               by="uids") %>% 
          				dplyr::mutate(model1=my_combn$result1[x], model2=my_combn$result2[x], combined = paste0(sign1," ", sign2)) %>% # add colums needed for later comparison
          				dplyr::count(model1, model2, combined) %>% 
          				tidyr::spread(key = combined, value = n)
      				}) %>% 
        			plyr::rbind.fill()
      		out[["sign"]] <- this_out
    	}
    	out_object <- get_make_results(object = object, 
                                   data = out, 
                                   metadata = NULL, 
                                   calc_type = rep("meta_network", length(out)),
                                   calc_info = rep(paste0("meta_network analysis with method",names(out)), each=length(out)),
                                   name = name) %>%
      				add_function_info(function_name = "meta_regression", 
                        params = list(result_index=result_index,
                                      method=method))
      	return(out_object)
	})