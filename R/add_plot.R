#' Generates a report from an metime_analyser object
#' @description Write an HTML or PDF report that summarizes one or more "me_time_analyser" objects and display all results
#' @param object a S4 object of class metime_analyser
#' @param result_item a plot or data table to be added as a result.
#' @param plot_data a dataframe corresponding to the result item.
#' @param calc_type a character defining the calc type of an object that has been added. Default set to NA.
#' @param calc_info a character defining the calc info of an object that has been added. Default set to NA.
#' @param result_index a character defining the result index of an object that has been added. Default set to 'result_index'
#' @param result_type  a character defining the result type of an object that has been added. Default set to 'result_type'
#' @param result_subtype  a character defining the result subtype of an object that has been added. Default set to 'result_subtype'
#' @param result_title  a character defining the result title of an object that has been added. Default set to 'result_title'
#' @return a S4 object of class metime_analyser
#' @export 
#
setGeneric("add_plot", function(object, result_item=NULL, plot_data=NULL, calc_type=NA, calc_info=NA, result_index = "result_index",result_type="result_type",result_subtype="result_subtype",result_title="result_title") standardGeneric("add_plot"))
setMethod("add_plot", "metime_analyser", function(object, result_item=NULL, plot_data=NULL, calc_type=NA, calc_info=NA, result_index = "result_index",result_type="result_type",result_subtype="result_subtype",result_title="result_title"){
  
  out <- object 

  # Error handling - check if parameter types are valid
  if(!all(sapply(c(calc_type,calc_info,result_index,result_type,result_subtype,result_title),length)==1)){
    warning("add_plot() calc_type, ,calc_info,result_index, result_type,result_subtype and result_title can only be a character"); run=F  
  }else{
    # Add result_item
    if(!is.null(out@results[[result_index]][["plots"]][[result_type]][[result_subtype]][[result_title]])){
      #Error handling - check if result_item already exists
      warning("add_results(): result_item already exists. Old data is overwritten.")
      result_nr <- which(names(out@results[[result_index]][["plots"]][[result_type]][[result_subtype]]) == result_title)
      out@results[[result_index]][["plots"]][[result_type]][[result_subtype]][[result_title]] <- result_item
      out@results[[result_index]][["information"]]$calc_type[result_nr] = calc_type
      out@results[[result_index]][["information"]]$calc_info[result_nr] = calc_info
    }else{
      out@results[[result_index]][["plots"]][[result_type]][[result_subtype]][[result_title]] <- result_item
      out@results[[result_index]][["information"]]$calc_type = c(out@results[[result_index]][["information"]]$calc_type,calc_type)
      out@results[[result_index]][["information"]]$calc_info = c(out@results[[result_index]][["information"]]$calc_info,calc_info)
    }
    
    # Add plot_data
    if(is.null(plot_data)){
      #no plot data availbale
    }else if(!is.null(out@results[[result_index]][["plots"]][[result_type]][[result_subtype]][[result_title]])){
      # Error handling - check if plot_data already exists
      warning("add_results(): plot_data already exists. Old data is overwritten.")
      object@results[[length(object@results)]]$plot_data[[result_index]] <- plot_data
    }else{
      result_nr <- which(names(object@results[[length(object@results)]]$plot_data) == result_index)
      object@results[[length(object@results)]]$plot_data[[result_index]] <- plot_data
      out@results[[result_index]][["information"]]$calc_type[result_nr] = calc_type
      out@results[[result_index]][["information"]]$calc_info[result_nr] = calc_info
    }
  }
  
  return(out)
  })

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
