#' Function to save interactive plots
#' @description extracts plot from plotter object and saves it as a widget
#' @examples save_plot_from_plotter(object)
#' @param object An object of class metime_plotter
#' @param out Character to specify path of the output file to save the widget in
#' @return saves the plot and returns nothing
#' @export
setGeneric("save_plot_from_plotter", function(object, out) standardGeneric("save_plot_from_plotter"))
setMethod("save_plot_from_plotter", "metime_plotter", function(object, out) {
			stopifnot(length(object@plot) == length(out))
			plots <- object@plot
			for(i in 1:length(out)) {
				saveWidget(plots[[i]], out[i])
			}
	})

