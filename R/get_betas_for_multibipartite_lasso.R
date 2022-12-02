#' Function to perform multibipartite style regression on a list of matrices
#' @description Performs multibipartite lasso in cv.glmnet style on a list of matrices that have
#' metabolite information from different platforms
#' @param list_of_mats a list with matrices and samples ordered similarly
#' @param alpha alpha for cv.glmnet regression. Defines style of penalty.
#' @param nfolds nfolds for cv.glmnet
#' @return returns a list with information of the combinations in context
#' @export
get_betas_for_multibipartite_lasso <- function(list_of_mats, # list of matrices that are divided based on platform or timepoints
					 alpha, # alpha parameter for glmnet
					 nfolds # nfolds parameter for glmnet
					 ) {
	#creating a list to store the data from glmnet
	#code exactly similar to the usual MLP 
	out <- list() # list to store the regression information for each metabolte
	count <- 1
	for(i in 1:length(list_of_mats)) {
		for(j in 1:length(list_of_mats)) {
			if(i != j) {
				x_mat <- as.matrix(list_of_mats[[i]])
				y_mat <- as.matrix(list_of_mats[[j]])
				fit_list <- list()
				for(k in 1:ncol(y_mat)) {
					fit_list[[k]] <- glmnet::cv.glmnet(x=x_mat, y=y_mat[,k], alpha=alpha, nfolds=nfolds)
					names(fit_list)[k] <- colnames(y_mat)[k]
				}
				out[[count]] <- fit_list
				names(out)[count] <- paste(names(list_of_mats)[j], names(list_of_mats)[i], sep="-")
				count <- count+1
			} else {
				fit_list <- list()
				mat <- as.matrix(list_of_mats[[i]])
				for(k in 1:ncol(mat)) {
					y <- as.matrix(mat[,k])
					x <- as.matrix(mat[,-k])
					fit_list[[k]] <- glmnet::cv.glmnet(x=x, y=y, alpha=alpha, nfolds=nfolds)
					names(fit_list)[k] <- colnames(mat)[k]
				} 
				out[[count]] <- fit_list
				names(out)[count] <- names(list_of_mats)[i]
				count <- count +1
			}

		}
	}
	return(out)
}

