#' Function to perform parametric tests on two sets of results or results of interest
#' @description Meta function to perform parametric tests on different types of results to perform meta analysis of interest
#' @param object An S4 object of class metime_analyser
#' @param results_index A character/numeric input for a results. The length should be equal to 1.
#' @param type character input to define the type of test. Two types are available: "t_test", "anova"
#' @param formula formula(class) input for the test functions. See rstatix::t_test() or rstatix::anova_test()
#' @param ... additional arguments for either rstatix::t_test() or rstatix::anova_test() based on input
#' @param name character input to set the name of the results
#' @seealso [meta_nonparametric_test], [rstatix::t_test], [rstatix::anova_test]
#' @returns An S4 object of class metime_analyser with this meta analysis stored as results
#' @export

setGeneric("meta_parametric_test", function(object, results_index, type, formula, name="meta_parametric_test_1", ...) standardGeneric("meta_parametric_test"))
setMethod("meta_parametric_test", "metime_analyser", function(object, results_index, type, formula, name="meta_parametric_test_1", ...) {
	if(!class(formula) %in% "formula") {
			warning("Formula should be of class formula. Exiting without making any changes.")
		}
		if(!type %in% c("t_test", "anova")) {
			warning("This type of test is not available. Exiting without making any changes.")
		}
		if(grep(name, names(object@results)) %>% length() >=1) {
            warning("name of the results was previously used, using a different name")
            index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
            index <- c(0:9)[grep(index, 0:9)+1]
            name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
        }
		results <- object@results[[results_index]]
		calc_type <- results$information$calc_type
		calc_info <- results$information$calc_info
		plot_data <- results$plot_data
		if(type %in% "anova") {
			plot_data <- lapply(seq_along(plot_data), function(ind) {
					test_res <- rstatix::anova_test(plot_data[[ind]], formula=formula, ...)
					return(test_res)
				})
		} else if(type %in% "t_test") {
			plot_data <- lapply(seq_along(plot_data), function(ind) {
					test_res <- rstatix::t_test(plot_data[[ind]], formula=formula, ...)
					return(test_res)
				})
		}
		exprs <- as.list(substitute(list(formula=formula, ...))[-1])
		str <- paste0(names(exprs), "=", sapply(exprs, function(x) {
		   if (is.character(x)) paste0('"', x, '"')
		   else as.character(x)
		 }), collapse = ", ")
		object <- object %>% 
				get_make_results(data=plot_data, metadata=NULL, calc_type=calc_type, calc_info=calc_info, name=name) %>%
				add_function_info(function_name="meta_parametric_test", 
					params=list(results_index=results_index, type=type, test_info=str))

})

