
#' Function to get the R environment
#' @description function to print the R environment
#' @return null
#' @export 
get_environment <- function() {
		x <- sapply(ls(), function(x) get(x, envir=.GlobalEnv))
}


