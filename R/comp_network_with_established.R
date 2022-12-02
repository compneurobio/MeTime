
#' Function to compare the network generated from data to an existing network
#' @include package_importer.R
#' @description function to perform fishers exact test to decide which network is the best.
#' @param calc_networks list of networks calculated from data
#' @param est_network established network to compare with the established network results
#' Make sure that this network has only two columns with names as node1 and node2
#' @return fisher test results with pval and test statistic
#' @export
setGeneric("comp_network_with_established", function(calc_networks, est_network) standardGeneric("comp_network_with_established"))
setMethod("comp_network_with_established", "metime_analyser", function(calc_networks, est_network) {
			for(i in 1:length(calc_networks)) {
				dummy <- calc_networks[[i]][,c("node1", "node2")]
				new_data <- generics::intersect(dummy, est_network)
				true_edges <- length(new_data[,1])


			}
	})



#Make RMarkdowns for networks with different kinds of tables and also fisher exact test
#Start making a write up with scientific questions
#check for courses on HELENA and choose a few
