#' register parallel backend
#' @description function to run in order to perform the analysis parallely thereby saving time
#' @param n_cores A number of specified cores.
#' @return set a parallel backend
#' @export


set_parallel_cores <- function(n_cores=NULL) {
  my_cores <- detectCores(all.tests = FALSE, logical = TRUE)
  if(!is.null(my_cores)) my_cores <- ifelse(n_cores>my_cores, my_cores, n_cores)
  else my_cores= my_cores - 1
  doParallel::registerDoParallel(cores=my_cores)
}
