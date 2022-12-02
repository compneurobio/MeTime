

#' Function to check if the data is already scaled or log transformed
#' @description Function to be applied on metime_analyser to check for log transformation and scaling
#' @param object An S4 object of metime_anlyser class
#' @param which_data the dataset/s to be checked
#' @returns NULL but checks if the data is scaled or not
#' @export
setGeneric("check_scaling_and_transformation", function(object, which_data) standardGeneric("check_scaling_and_transformation")) 
setMethod("check_scaling_and_transformation", "metime_analyser", function(object, which_data) {
			list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
			test <- lapply(names(list_of_data), function(x) {
					data <- list_of_data[[x]]
					y <- apply(data, 2, mean)
					z <- apply(data, 2, sd)
					y <- y==0
					z <- z==1
					a <- apply(data, 2, function(m) {
							m <- m > 0
							return(m)
						})
					if(all(y) && all(z)) {
						print(paste(x, "dataset is scaled", sep=" "))
					} else {
						warning(paste(x, "dataset is not scaled", sep=" "))
					}
					if(all(a)) {
						warning(paste(x, "dataset is neither scaled nor log transformed", sep=" "))
					} else if(!all(a) && !all(y) && !all(z)) {
						warning(paste(x, "dataset is NOT scaled but could be log transformed", sep=" "))
					}
				})
	})
