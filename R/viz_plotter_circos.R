

#' Setting up standard wrapper for all circos plots for a metime_plotter object. 
#' @description plot function for metime_plotter object with different inputs to specialize plots. Used for all calc outputs.
#' @param object S4 object of class metime_plotter
#' @param aesthetics list for aesthetics. eg: list(list(x="colname",y="colname",color="colname", shape="colname"), list(...)) for "dot" plot and "heatmap"
#' plot, for heatmap: list(x="colname", y="colname", fill="colname"). Additionally two other character vectors are allowed namely .$vis and .$strats for text
#' and for facet wrapping.
#' @export
setGeneric("viz_plotter_circos", function(object, aesthetics, outfile, layout_by) standardGeneric("viz_plotter_circos")) 
setMethod("viz_plotter_circos", "metime_plotter", function(object, aesthetics, outfile, layout_by) {
			pdffile <- pdf(outfile)

			dev.off()
	}) 

