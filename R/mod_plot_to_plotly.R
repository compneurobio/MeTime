#' Function to change ggplot object to ggplotly object with text
#' @description converts ggplot to plotly 
#' @param object An S4 object of class metime_analyser
#' @param results_index index of results
#' @returns metime_analyser object with updated plot
#' @export
setGeneric("mod_plot_to_plotly", function(object, results_index) standardGeneric("mod_plot_to_plotly"))
setMethod("mod_plot_to_plotly", "metime_analyser", function(object, results_index) {
		
	})