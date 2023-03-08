#' Function that adds a result to a metime anal
#' @description Function to add information about the method applied to the dataset
#' @param object S4 object of class metime_analyser
#' @param x a data.frame or plot that should be added to the list of results
#' @param functions_applied a character providing information on the calculation
#' @param calc_info a character that provides information on the calculation
#' @param calc_type a character that defines the type of calculation
#' @param results_index character or numeric of length=1 to define the results to which you want to add x.
#' Default is set to NULL and that implies that new results are generated 
#' @return Saves the report as html/pdf
#' @export 

setGeneric("add_result", function(object, x, name=NULL, functions_applied=NULL, calc_info=NULL, calc_type=NULL, results_index=NULL) standardGeneric("add_result"))
setMethod("add_result", "metime_analyser", function(object, x, name=NULL, functions_applied=NULL, calc_info=NULL, calc_type=NULL, results_index=NULL){
  if(is.null(name)) name <- paste0("result_",length(object@results)+1)
  if(is.null(results_index)) {
      object@results[[name]] <- list(
        functions_applied = list(
          add_result = ifelse(is.null(functions_applied), "no information available",as.character(functions_applied))
        ),
        plot_data = list(data.frame()),
        information =list(calc_type=ifelse(is.null(calc_type), "no information available",as.character(calc_type)),
                      calc_info=ifelse(is.null(calc_info), "no information available",as.character(calc_info))
                      ),
        plots = list(
            list(
              name = list(
                name=x
              )
            )
        )
      )
      names(object@results[[name]]$plots[[1]]) <- name
      names(object@results[[name]]$plots[[1]][[name]])<- name
      names(object@results[[name]]$functions_applied$add_result) <- paste0("par",1:length(object@results[[name]]$functions_applied$add_result)) 
  } else {
      results <- object@results[[results_index]]
      results$plots[[length(results$plots)+1]] <- list(name=list(name=x))
      results$functions_applied[[length(results$functions_applied)+1]] <- c(addition="new_plot") 
      names(results$functions_applied)[length(results$functions_applied)] <- "add_result"
      object@results[[results_index]] <- results
  }
  
  return(object)
})
