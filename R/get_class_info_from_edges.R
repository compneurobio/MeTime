#' Get summary on class edges 
#' @description Function to check how the different edges in a GGM are associated to their respective classes(it could be super-pathway or sub-pathway)
#' @param calc_networks a list of calculated networks
#' @param metadata a dataframe containing the metadata of the edges present.
#' @param phenotypes a character vector to define phenotypes that were used for correcting the data
#' @return A dataframe with information on different type of edges present
#' @export
get_class_info_from_edges <- function(calc_networks, metadata, phenotypes) {
					out <- list()
					for(i in 1:length(calc_networks)) {
								rm_phen <- phenotypes
								calc_networks[[i]] <- calc_networks[[i]][!calc_networks[[i]]$node1 %in% rm_phen, ]
								calc_networks[[i]] <- calc_networks[[i]][!calc_networks[[i]]$node2 %in% rm_phen, ]
								network <- calc_networks[[i]]
								network$node1 <- metadata[as.character(calc_networks[[i]]$node1) %in% metadata$name, "group"]	
								network$node2 <- metadata[as.character(calc_networks[[i]]$node2) %in% metadata$name, "group"]
								
								network <- network[ ,c("node1", "node2")]
								network <- data.table::setDT(network)[ ,list(count=.N), names(network)]
								network <- na.omit(network)
								out[[i]] <- network
					}
					return(out)
	}


