#' Function that generates a report from an metime_analyser object
#' @description Function to add information about the method applied to the dataset
#' @param object S4 object of class metime_analyser
#' @param title a character string used as title of the report. By default is set to metime + system time.
#' @param file a character string used as file name of the output report and table. By default is set to metime + system time.
#' @param table a logical indicating whether you want to save the plot data as an xlsx file. If TRUE an xlsx file with the file name specified by the title will be saved.
#' @param device a character string specifying the format of the report. By default is set to html. Other options include pdf.
#' @param interactive a logical indicating whether plots are interactive - only possible for html. Default set to FALSE
#' @param author a character specifying the name of the author
#' @return Saves the report as html/pdf
#' @export 

setGeneric("get_report", function(object, title=NULL, file=NULL, table=F, device="html", interactive=F, author=NULL) standardGeneric("get_report"))
setMethod("get_report", "metime_analyser", function(object, title=NULL, file=NULL, table=F,device="html", interactive=F, author=NULL){
  results <- object@results
  # sanity checks ####
  get_header <- function(title, type="html"){
    # function that generates the heading of the rmarkdown document
    if(type == "html"){
      # make HTML report
      out <- c(paste0('---\n',
                      'title: "', title,'"\n',
                      'date: ',Sys.Date(),'\n',
                      ifelse(!is.null(author), paste0('author: ',paste0(author, collapse=", "),'\n'),''),
                      'output:\n',
                      " rmdformats::readthedown:\n  highlight: kate\n  code_folding: hide\n--- \n"))
    }else if(type == "pdf"){
      out <- paste0('---\n',
                      'title: "', title,'"\n',
                      'date: ',Sys.Date(),'\n',
                      ifelse(!is.null(author), paste0('author: ',paste0(author, collapse=", "),'\n'),''),
                      'output: \npdf_document:\n',
                          'toc: true\n',
                           'toc_depth: 2\n',
                    '---\n')
    }else{
      out=NULL
    }
    return(out)
  }
  
  get_analyzer_info <- function(object){
    # function that prints all information of the analyzer object in rmarkdown
    out <- paste0("# Analyzer info\n```{r,echo=F, warning=FALSE, message=FALSE}\n",
             'results <- object@results\n ',
             'print(object)',
             "\n```\n", collapse="")
    return(out)
  }
  
  # additional functions ####
  get_p<-function(...){
    # function that makes a new paragraph section in rmarkdown
    p_list <- list(...)
    out <- lapply(p_list,function(x)paste0('\n',x,'\n'))
    out <- unlist(out)
    return(out)
  }
  
  get_heading <- function(text, level=1, tabset=F, type="html"){
    if(type=="html"){
      out <- paste0('\n', paste0(rep("#",each=level), collapse="")," ", text, ifelse(tabset, "{.tabset}",""),'\n')
    }
    else if(type=="pdf"){
      out <- paste0('\n', paste0(rep("#",each=level), collapse="")," ", text,ifelse(level==1,paste0(" {#",text,"}"),""), '\n')
    }
    else{
      out <- ""
    }
    return(out)
  }
  
  get_table <- function(table_location, caption=NULL, rownames=F, type="html"){
    if(type=="html"){
      out <- paste0('DT::datatable(',table_location,', caption = "',caption,'", rownames=',rownames, ',' ,
      'extensions = "Buttons,\n options = list(dom = "Blfrtip,\nbuttons = c("csv", "excel")))')
    }
    else if(type=="pdf"){
      out <- paste0('knitr::kable(',table_location,'caption = "',caption,'", rownames=',rownames,'")')
    }
    else{
      out <- ""
    }
    return(out)
  }
  
  get_plot <- function(plot_location, title=NULL, caption=NULL, type="html", interactive=F){
    if(type=="html" & interactive){
      out <- c(ifelse(!is.null(title), paste0('\n',title, ' \n'), ""), # add title
               "```{r,echo=F, warning=FALSE, message=FALSE}\n",
               "plotly::ggplotly(",
               plot_location,
               ")",
               "\n```",
               ifelse(!is.null(caption), paste0('\n',caption, ' \n'), "") # add caption
      )
    }
    else if(type=="html" & !interactive){
      out <- c(ifelse(!is.null(title), paste0('\n',title, ' \n'), ""), # add title
               "```{r,echo=F, warning=FALSE, message=FALSE}\n",
               plot_location,
               "\n```",
               ifelse(!is.null(caption), paste0('\n',caption, ' \n'), "") # add caption
      )
    }
    else if(type=="pdf"){
      out <- c(ifelse(!is.null(title), paste0('\n',title, ' \n'), ""), # add title
               "```{r,echo=F, warning=FALSE, message=FALSE}\n",
               plot_location,
               "\n```",
               ifelse(!is.null(caption), paste0('\n',caption, ' \n'), "") # add caption
      )
    }
    else{
      out <- ""
    }
    return(out)
  }
  
  get_fun_info <- function(fun, type="html"){
    get_package <- utils::find(fun) %>% gsub(pattern="package:",replacement = "")
    if(length(get_package)==0){
      out <- "NA"
    }
    else if(get_package == ".GlobalEnv"){
      out <- "NA"
    }
    else{
      get_db <- tools::Rd_db(package = get_package)
      if(length(grep(pattern=fun, names(get_db), value=T))>0){
        my_fun_db <- get_db[[grep(pattern=fun, names(get_db), value=T)]]
        out <- capture.output(tools::Rd2HTML(my_fun_db))[-c(1:7)]
      }
      else out <- "function not found"
      
      out <- out %>% 
        gsub(pattern=fun,fixed=T,replacement = "
         ") 
    }
    out <- out %>% gsub(pattern='"', replacement = "'", fixed = T)
    return(out)
  }
  
  get_pipe <- function(list, prefix=NULL,type="html"){
    if(type=="html"){
      #gather pipe information
      pipe_info <- lapply(seq_along(list), function(x){
        this_info <- data.frame(var = names(list[[x]]) %>% gsub(pattern="[0-9]", replacement = ""), par=as.character(list[[x]]), stringsAsFactors = F)
        my_params <-paste0("(",
                           lapply(unique(this_info$var), function(y){
                             if(y!="mutations"){
                               # if mutations then this item was used in mod_mutate or mod_filter
                               paste0(y, " = ", ifelse(length(which(this_info$var==y))==1,
                                                       paste0("'",this_info$par[which(this_info$var==y)],"'"), 
                                                       paste0("c('",paste0(this_info$par[which(this_info$var==y)], collapse="', '"),"')")))
                             }else{
                               paste0( ifelse(length(which(this_info$var==y))==1,
                                              paste0(this_info$par[which(this_info$var==y)],"'"), 
                                              paste0("c('",paste0(this_info$par[which(this_info$var==y)], collapse=" , "),")")))
                             }
                             })%>%
            paste0(collapse = ", "), ")") %>% gsub(pattern='"',replacement = "'",fixed=T)
        paste0('   <button class=button_fun data-toggle=modal data-target=',paste0("#",prefix,x),' >', names(list)[x],'</button><span>',my_params,' %>% </span><br>')
      })
      pipe_info[[length(pipe_info)]] <- pipe_info[[length(pipe_info)]] %>% gsub(pattern="%>%", replacement = "", fixed=T)
      
      pipe_modal <- lapply(seq_along(list), function(x){
        paste0("```{r ,echo=F}\n",
               'bsplus::bs_modal(id = "',paste0(prefix,x),'",',
               'title = "',names(list)[x],'",',
               'size = "large",',
               'body = htmltools::HTML("',paste0(get_fun_info(names(list)[x]),collapse=""),'")',
               ')\n``` \n')
      }) %>% unlist()
      
      out <- c(paste0("\n```{r, echo=F}\n", 
                    'htmltools::HTML("<div class=body_fun> object %>% <br>',
                    '<div style= margin-left:5px>',
                    paste0(pipe_info, collapse=""), 
                    '</div>',
                            '</div>")\n```\n'
                    ),
               pipe_modal
               )
    }
    else if(type=="pdf"){
      out <- lapply(seq_along(list), function(x){
        this_info <- data.frame(var = names(list[[x]]) %>% gsub(pattern="[0-9]", replacement = ""), par=as.character(list[[x]]), stringsAsFactors = F)
        my_params <-paste0("(",lapply(unique(this_info$var), function(y) 
          paste0(y, " = ", ifelse(length(which(this_info$var==y))==1,
                                paste0("'",this_info$par[which(this_info$var==y)],"'"), 
                                paste0("c('",paste0(this_info$par[which(this_info$var==y)], collapse="', '"),"')")))) %>%
            paste0(collapse = ", "), ")") %>% gsub(pattern='"',replacement = "'",fixed=T)
        out <- paste0("\n >> ###### ", names(list)[x],"", my_params," *%>%*\n")
      })
      out[[length(out)]] <- out[[length(out)]] %>% gsub(pattern="*%>%*", replacement = "", fixed=T)
      out <- paste0("\n > ###### object *%>%* \n", paste0(out, collapse=""))
    }else{
      out = ""
    }
    return(out)
  }
  
  # specify outputs ####
  out_file <- ifelse(is.null(file), paste0("metime",Sys.Date()), tools::file_path_sans_ext(file))
  out_title <- ifelse(is.null(title),paste0("metime",Sys.Date()),tools::file_path_sans_ext(file))
  
  
  #### header ####
  out_header <- get_header(title=ifelse(is.null(title),paste0("metime",Sys.Date()),tools::file_path_sans_ext(file)),
                           type=device)
  out_analyzer <- get_analyzer_info(object = object)
  # CSS block 
  out_css <- ifelse(device=="html", 
  "```{css, echo=FALSE}
#content{max-width:100%}
#sidebar{bacgkround:black}
#postamble{border-top:solid 10px black}
.button_fun{padding: 0 2px; border:transparent; background:transparent; color:#9F2042; }
.button_fun:hover{padding: 0 2px; border:transparent; background:transparent; color:#9F2042;font-weight:bold}
.body_fun{border-radius: 10px; background:rgb(214, 214, 214); padding:10px; margin:10px 0; }
```",
  "")
  # JS block
  out_js <- ''
  
  #### body ####
  table_list <- list()
  table_info <- data.frame()
  out_body <- lapply(names(results), function(section_name){
    # add section title
    section_title = get_heading(text=section_name, tabset=T, level=1, type=device)
    section_sub_title = get_p(#'\n',"General information",'\n',
                              #'\n',paste0("Calculation type:  ", results[[section_name]][["information"]]$calc_type)%>% gsub(pattern='"', replacement = "'", fixed=T),
                              '\n',paste0("Information:  ", results[[section_name]][["information"]]$calc_info)%>% gsub(pattern='"', replacement = "'", fixed=T),
                              '\n', '\n')
    section_pipe = get_pipe(results[[section_name]][["functions_applied"]], prefix=section_name, type=device)
    
    # add tabs to section
    section_tab = lapply(seq_along(results[[section_name]]$plots), function(tab_nr){
      # plot all subitems
      lapply(seq_along(results[[section_name]][["plots"]][[tab_nr]]), function(section_plots){
        if(all(class(results[[section_name]][["plots"]][[tab_nr]][[section_plots]])=="list")){
          lapply(seq_along(results[[section_name]][["plots"]][[tab_nr]][[section_plots]]), function(section_sub_plots){
            if(is.data.frame(results[[section_name]][["plots"]][[tab_nr]][[section_plots]][[section_sub_plots]])){
              c(get_heading(text=names(results[[section_name]][["plots"]][[tab_nr]])[section_plots], level=2, tabset=F, type=device), # add plot name
                get_table(plot_location = paste0('results[["',section_name,'"]][["plots"]][[',tab_nr,']][[',section_plots,']][[',section_sub_plots,']]'), type=device)
              )
            }else{
              #add a plot to the report
              c(get_heading(text=names(results[[section_name]][["plots"]][[tab_nr]])[section_plots], level=2, tabset=F, type=device), # add plot name
                get_plot(plot_location = paste0('results[["',section_name,'"]][["plots"]][[',tab_nr,']][[',section_plots,']][[',section_sub_plots,']]'), type=device, interactive = interactive)
              )
            }
          }) %>% unlist()
        }else{
          c(get_heading(names(results[[section_name]][["plots"]][[tab_nr]])[section_plots], level=2, tabset=F, type=device), # add plot name
            get_plot(plot_location = paste0('results[["',section_name,'"]][["plots"]][[',tab_nr,']][[',section_plots,']]'), type=device)
          )
        }
        # add all plots
      }) %>% unlist()
    }) %>% unlist()
    # save_plot_data 
    if(table){
      table_list<- append(table_list, results[[section_name]][["plot_data"]])
      names(table_list) = NULL
      table_info = rbind(table_info, data.frame(calc_type=results[[section_name]]$information$calc_type,
                                                calc_info=results[[section_name]]$information$calc_info))
    }
    
    c(
      section_title,
      section_sub_title,
      section_pipe,
      #section_info,
      section_tab
    )
  }) %>% unlist()
  
  # Footer ####
  out_footer <- "# Session info\n```{r session_info,include=T, warning=FALSE, message=FALSE}\nprint(sessionInfo())\n```\n"
  
  # Write xlsx
  if(table){
    name(table_list)
    my_info$sheet_nr <- 1:nrow(my_info)
    xlsx::write.xlsx(x=my_info, file = paste0(out_file,".xlsx"),sheetName = "info")
    
    for(x in seq_along(my_tables)){
      xlsx::write.xlsx(x=my_tables[[x]], file = paste0(out_file,".xlsx"),sheetName = paste0("sheet",my_info$sheet_nr[x]), append = T)
    }
  }
  # write report ####
  writeLines(
    text=c(
      out_header,
      out_css,
      out_js, 
      out_analyzer,
      out_body,
      out_footer), # html code
    paste0(out_file,".rmd") #path for saving the file
  )
  
  # convert rmd to html
  output_format = NULL
  if(device %in% c("pdf","word")) output_format <- paste0(device, "_document")
  rmarkdown::render(input=paste0(out_file,".rmd"),output_format = output_format)
})
