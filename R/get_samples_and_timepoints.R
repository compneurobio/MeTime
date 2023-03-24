#' Summarize number of time points and the total number of samples available at that point
#' @description Obtain the number of unique samples available at each time point.
#' @param objecta a S4 object of the class "metime_analyzer".
#' @param which_data a character to define which dataset is to be used.
#' @return A dataframe with time points and number of samples at each time point.
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


