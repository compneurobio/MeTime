#' Function to save plot_data from an analyser_objec into an xlsx file
#' @description Saves plot_data from an analyser_objec into an xlsx file
#' @param object An object of class metime_analyser
#' @param file a character defining the path where to save the xlsx file
#' @return saves all plot data in an xlsx sheet
#' @export
setGeneric("write_results", function(object, file=NULL) standardGeneric("write_results"))
setMethod("write_results", "metime_analyser", function(object, file=NULL) {
  
  out_file <- ifelse(is.null(file), paste0("metime",Sys.Date()), tools::file_path_sans_ext(file))

  out_tables <- lapply(seq_along(object@results), function(result_nr){
    ##iterate across results items
    lapply(seq_along(object@results[[result_nr]]$plot_data),function(data_nr){
      ## iterate within results items
      pipe_info <- lapply(seq_along(object@results[[result_nr]]$functions_applied), function(x){
        ## get pipe information
        this_info <- data.frame(var = names(object@results[[result_nr]]$functions_applied[[x]]) %>% gsub(pattern="[0-9]", replacement = ""), par=as.character(object@results[[result_nr]]$functions_applied[[x]]), stringsAsFactors = F)
        my_params <-paste0("(",lapply(unique(this_info$var), function(y) 
          paste0(y, " = ", ifelse(length(which(this_info$var==y))==1,
                                  paste0("'",this_info$par[which(this_info$var==y)],"'"), 
                                  paste0("c('",paste0(this_info$par[which(this_info$var==y)], collapse="', '"),"')")))) %>%
            paste0(collapse = ", "), ")") %>% gsub(pattern='"',replacement = "'",fixed=T)
        paste0(names(object@results[[result_nr]]$functions_applied)[x],
               "", 
               my_params,
               ifelse(x == length(seq_along(object@results[[result_nr]]$functions_applied)),""," %>% ")
        )
      }) %>% unlist() %>% paste0(collapse="")
      
      list(
        plot_data = object@results[[result_nr]]$plot_data[[data_nr]],
        result_name = names(object@results[[result_nr]]$plot_data)[data_nr],
        calc_type = object@results[[result_nr]]$information$calc_type[data_nr],
        calc_info = object@results[[result_nr]]$information$calc_info[data_nr],
        pipe_info = pipe_info
      )
    })
  }) %>% Reduce(c,.)
  
  info_table <- data.frame(
    result_name = sapply(out_tables,function(x) x[["result_name"]]) %>% as.character(),
    sheet_name = paste0("result_",1:length(out_tables)),
    calc_type = sapply(out_tables,function(x) x[["calc_type"]]) %>% as.character(),
    calc_info = sapply(out_tables,function(x) x[["calc_info"]]) %>% as.character(),
    pipe_info = sapply(out_tables,function(x) x[["pipe_info"]]) %>% as.character(),
    stringsAsFactors = F
  )
  
  # write xslx file with info_table
  xlsx::write.xlsx(info_table, 
                   file=paste0(out_file, ".xlsx"), 
                   sheetName="Information", 
                   row.names=FALSE)
  ## add all other sheets
  for(sheet_nr in seq_along(out_tables)){
    xlsx::write.xlsx(out_tables[[sheet_nr]]$plot_data, 
                     file=paste0(out_file, ".xlsx"), 
                     sheetName=info_table$sheet_name[sheet_nr], 
                     append = T,
                     row.names=FALSE)
  }
	}) 
