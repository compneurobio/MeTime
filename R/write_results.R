#' Function to save plot_data from an analyser_objec into an xlsx file
#' @description Saves plot_data from an analyser_objec into an xlsx file (no Java)
#' @param object An object of class metime_analyser
#' @param file a character defining the path where to save the xlsx file
#' @return saves all plot data in an xlsx sheet
#' @seealso [write_report], [write_data]
#' @export
setGeneric("write_results", function(object, file=NULL) standardGeneric("write_results"))
setMethod("write_results", "metime_analyser", function(object, file=NULL) {

  if(!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Please install 'openxlsx' to use this write_results implementation.")
  }

  out_file <- ifelse(is.null(file), paste0("metime", Sys.Date()), tools::file_path_sans_ext(file))

  out_tables <- lapply(seq_along(object@results), function(result_nr){
    ## iterate across results items
    lapply(seq_along(object@results[[result_nr]]$plot_data), function(data_nr){
      ## iterate within results items
      list(
        plot_data = object@results[[result_nr]]$plot_data[[data_nr]],
        result_name = names(object@results[[result_nr]]$plot_data)[data_nr],
        calc_type = object@results[[result_nr]]$information$calc_type[data_nr],
        calc_info = object@results[[result_nr]]$information$calc_info[data_nr]
      )
    })
  }) %>% Reduce(c, .)

  info_table <- data.frame(
    result_name = sapply(out_tables, function(x) x[["result_name"]]) %>% as.character(),
    sheet_name = paste0("result_", 1:length(out_tables)),
    calc_type = sapply(out_tables, function(x) x[["calc_type"]]) %>% as.character(),
    calc_info = sapply(out_tables, function(x) x[["calc_info"]]) %>% as.character(),
    stringsAsFactors = FALSE
  )

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Information")
  openxlsx::writeData(wb, "Information", info_table)

  for(sheet_nr in seq_along(out_tables)) {
    sheet_name <- info_table$sheet_name[sheet_nr]
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, out_tables[[sheet_nr]]$plot_data)
  }

  openxlsx::saveWorkbook(wb, file = paste0(out_file, ".xlsx"), overwrite = TRUE)
})
