#' Function to map edge lists to HMDB and KEGG IDs
#' @description Function to map metabolite names 
#' @param object An S4 object of class metime_analyser
#' @param results_index Name/index for the network results
#' @param metanno Dataframe with mappings
#' @returns mapped results
#' @export
setGeneric("meta_fisher_for_networks", function(object, results_index, metanno, app_lookup) standardGeneric("meta_fisher_for_networks"))
setMethod("meta_fisher_for_networks", "metime_analyser", function(object, results_index, metanno, app_lookup) {
		results <- object@results[[results_index]]
		data <- results$plot_data$edge
		data <- tidyr::pivot_wider(data, names_from=node1, values_from=values) %>% as.data.frame()
		#app_lookup <- feather::read_feather("app_lookup.feather")

		# Extract relevant RECON ids ----

		# get recon ids for metabolites using HMDB
		recon_ids_HMDB <- app_lookup %>% 
  			dplyr::inner_join(metanno, by = c("id" = "HMDB"))
		# remove duplicated entries
		recon_ids_HMDB <- recon_ids_HMDB[!duplicated(recon_ids_HMDB$recon_id) & !duplicated(recon_ids_HMDB$id) & !duplicated(recon_ids_HMDB$BIOCHEMICAL),]

		# get recon ids for metabolites using KEGG
		recon_ids_KEGG <- app_lookup %>% 
  			dplyr::inner_join(metanno, by = c("id" = "KEGG"))
		# remove duplicated entries
		recon_ids_KEGG <- recon_ids_KEGG[!duplicated(recon_ids_KEGG$recon_id) & !duplicated(recon_ids_KEGG$id) & !duplicated(recon_ids_KEGG$BIOCHEMICAL),]
		# remove entries covered by HMDB already
		recon_ids_KEGG %<>%
  			dplyr::filter(!(recon_id %in% recon_ids_HMDB$recon_id))

		# merge the two
		recon_ids <- recon_ids_HMDB %>%
  			dplyr::full_join(recon_ids_KEGG)

		# Get Distance Matrix ----

		# subset distance matrix to mapped recon ids
		met_ids <- colnames(out_distances) %>%
  			intersect(recon_ids$recon_id) 

		# cut distance matrix to only genes and metabolites of interest
		distmat <- out_distances[met_ids, met_ids]  

		# reassign metabolite names instead of recon ids
		names <- data.frame(recon_id=colnames(distmat)) %>%
  			dplyr::left_join(recon_ids[,c("recon_id","BIOCHEMICAL")]) 
		names <- names[!duplicated(names$recon_id),]

		# use metabolite names as column and row names
		colnames(distmat) <- rownames(distmat) <- names$BIOCHEMICAL
		return(distmat)
	})