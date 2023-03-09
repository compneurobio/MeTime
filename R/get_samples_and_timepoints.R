
#' Function to know the number of timepoints and the total number of samples available at that point
#' @description A method applied onto s4 object of class "metime_analyser" so as to obtain the number of unique samples available
#' at each timepoint. 
#' @examples
#' # newdata <- get_samples_and_timepoints(object=metime_analyser_object, which_data="Name of dataset of interest")
#' @param object An object of class metime_analyser
#' @param which_data Name of the dataset in context
#' @return A data table with timepoints and number of samples at each timepoint
#' @export
setGeneric("get_samples_and_timepoints", function(object, which_data) standardGeneric("get_samples_and_timepoints"))

setMethod("get_samples_and_timepoints", "metime_analyser", function(object, which_data){
		data <- object@list_of_row_data[names(object@list_of_data) %in% which_data]
		data <- data[[1]]
		unique_timepoints <- rownames(data) %>% gsub(pattern="[a-z|A-Z][0-9]+_", replacement="") %>% unique()
		levels <- unique_timepoints %>% gsub(pattern="[a-z|A-z]", replacement="") %>% as.numeric() %>% order()
		unique_timepoints <- unique_timepoints[levels]
		samples_count <- vapply(unique_timepoints, function(x) {
						index <- x %in% data$time 
						index <- index[index==TRUE]
						return(length(index))
					}, numeric(1))
		out <- as.data.frame(cbind(as.character(unique_timepoints), samples_count))
		colnames(out) <- c("timepoints", "number of samples")
		return(out)
	})


