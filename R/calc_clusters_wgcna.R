#' Calculate cluster assignment from WGCNA
#' @description calculate cluster assignment to specified col_data from Weighted Correlation Network Analysis
#' @param object a S4 object of the class metime_analyzer.
#' @param which_data a character to define which dataset is to be used.
#' @param baseline a character to define the baseline time point which is used for cluster calculation.
#' @param cols_for_meta A list of named character vector of length equal to which_data. 
#' Example: list(dataset1=c(new_name1="colname", new_name2="colname"))
#' Default is set to NULL
#' @param name Name of the results. Default is set to "WGCNA_clusters_1"
#' @param soft_power Numeric soft-thresholding power to use. When `NULL`, the
#' first candidate power with a truncated R-squared greater than 0.8 is used.
#' If no candidate meets that threshold, the function warns and returns without
#' changing `object`. Supply a numeric value to override the automatic choice.
#' @param ... multiple parameters separated by commas passed to cutreeDynamic (R package: dynamicTreeCut). 
#' Parameters include: minClusterSize, pamDendroRespect ...
#' @details Based on the method described in the WGCNA tutorials for step-by-step network construction and module detection. The idea is to find
#' modules of metabolites at baseline. This method can be found in detail here [WGCNA: an R package for weighted correlation network analysis](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)
#' @seealso [dynamicTreeCut::cutreeDynamic], [get_metadata_for_columns], [mod_trans_eigendata]
#' @return a S4 object of class "metime_analyser" with cluster information appended to col_data of which_data
#' @export
setGeneric("calc_clusters_wgcna", function(object, which_data, baseline="t0", cols_for_meta=NULL, name="WGCNA_clusters_1", soft_power=NULL, ...) standardGeneric("calc_clusters_wgcna"))
setMethod("calc_clusters_wgcna", "metime_analyser", function(object, which_data, baseline="t0", cols_for_meta=NULL, name="WGCNA_clusters_1", soft_power=NULL, ...) {
        # Sanity checks
        if(length(which_data)!=1 & !which_data %in% names(object@list_of_data)) warning("calc_clusters_wgcna(): which_data not in metime_analyzer, or more than one which_data selected")
        else if(length(grep(baseline, rownames(object@list_of_data[[which_data]])))<=0) warning("calc_clusters_wgcna(): baseline timepoint not found")
        else{
                while(grep(name, names(object@results)) %>% length() >=1) {
                        warning("name of the results was previously used, using a different name")
                        index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
                        index <- c(0:9)[grep(index, 0:9)+1]
                        name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
                }
                # Get data
                data <- object@list_of_data[[which_data]]
                col_data <- object@list_of_col_data[[which_data]]
                data <- data[grep(baseline, rownames(data)), ]
                lnames <- names(data)
                # Choose a set of soft-thresholding powers
                powers = c(c(1:10), seq(from = 12, to=20, by=2))
                # Call the network topology analysis function
                sft = WGCNA::pickSoftThreshold(data, powerVector = powers) # verbose 0 = silent, verbose 5 = print info
                if(is.null(soft_power)) {
                        eligible_powers <- sft$fitIndices$Power[
                                which(sft$fitIndices$truncated.R.sq > 0.8)
                        ]
                        if(length(eligible_powers) == 0) {
                                warning(
                                        paste0(
                                                "calc_clusters_wgcna(): no candidate soft power has ",
                                                "truncated R-squared > 0.8; stopping without changes. ",
                                                "Provide soft_power as a numeric value if you wish to override this selection."
                                        )
                                )
                                return(object)
                        }
                        soft_power <- eligible_powers[[1]]
                        message(
                                "Using soft_power = ", soft_power,
                                ", the first candidate with truncated R-squared > 0.8. ",
                                "Provide soft_power as a numeric value to choose a different power."
                        )
                } else if(!is.numeric(soft_power) || length(soft_power) != 1 ||
                          !is.finite(soft_power) || soft_power <= 0) {
                        stop("calc_clusters_wgcna(): soft_power must be NULL or one positive, finite numeric value.")
                }
                adjacency <- WGCNA::adjacency(data, power=soft_power)
                # Turn adjacency into topological overlap
                TOM = WGCNA::TOMsimilarity(adjacency)
                dissTOM = 1-TOM
                # Call the hierarchical clustering function
                geneTree = hclust(as.dist(dissTOM), method = "average")  
                dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM, ...)
                modules <- cbind.data.frame(lnames, dynamicMods)
                modules <- modules[order(modules$lnames), ]
                colnames(modules) <- c("id", "modules")
                rownames(modules) <- modules$id
                if(is.null(cols_for_meta)) {
                        metadata <- NULL
                } else {
                        metadata <- get_metadata_for_columns(object=object, 
                                which_data=which_data, columns=cols_for_meta)
                }
                object <- get_make_results(object=object, data=list(modules=modules), metadata=metadata, calc_type="clusters", 
                        calc_info=paste("WGCNA clusters info of", which_data, sep=""), name=name)
                object <- add_function_info(object=object, 
                                            function_name="calc_clusters_wgcna", 
                                            params=list(which_data=which_data, baseline=baseline, name=name, soft_power=soft_power, ...))
        }
        return(object)
})
