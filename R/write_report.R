#' Generates a report from an metime_analyser object
#' @description Write an HTML or PDF report that summarizes one or more "me_time_analyser" objects and display all results
#' @param object a S4 object of class metime_analyser
#' @param title a character string used as title of the report. By default is set to metime + system time.
#' @param file a character string used as file name of the output report and table. By default is set to metime + system time.
#' @param write_results a logical indicating whether you want to save the plot data as an xlsx file. If TRUE an xlsx file with the file name specified by the title will be saved.
#' @param device a character string specifying the format of the report. By default is set to html. Other options include pdf.
#' @param interactive a logical indicating whether plots are interactive - only possible for html. Default set to FALSE
#' @param author a character specifying the name of the author
#' @return Saves the report as html/pdf
#' @export 

setGeneric("write_report", function(object, title=NULL, file=NULL, write_results=F, device="html", interactive=F, author=NULL) standardGeneric("write_report"))
setMethod("write_report", "metime_analyser", function(object, title=NULL, file=NULL, write_results=F,device="html", interactive=F, author=NULL){
<<<<<<< Updated upstream
  
  # set parameters -------
  ## set name of the file
  out_file <- ifelse(is.null(file), paste0("metime",Sys.Date()), tools::file_path_sans_ext(file))
  
  # process data -------
  ## save analyzer objects as a list
  if(!is.list(object)){ 
    analyzer_list <- list(object)
    analyzer_list_names <- "object1"
  }else{
    analyzer_list <- object
  }
  
  ## reorder the results list and annotate as results_report
  results_list <- list()
  for(nr_analyzer in seq_along(analyzer_list)){
    results_list[[nr_analyzer]] <- lapply(seq_along(analyzer_list[[nr_analyzer]]@results),function(nr_results){
      list(
        functions_applied = analyzer_list[[nr_analyzer]]@results[[nr_results]]$functions_applied,
        items = lapply(seq_along(analyzer_list[[nr_analyzer]]@results[[nr_results]]$plot_data),function(nr_data){
          list(
            h2=names(analyzer_list[[nr_analyzer]]@results[[nr_results]]$plot_data)[nr_data], # name of item
            plot_data = analyzer_list[[nr_analyzer]]@results[[nr_results]]$plot_data[[nr_data]], #plot data of item
            calc_type = analyzer_list[[nr_analyzer]]@results[[nr_results]]$information$calc_type[nr_data], #calc type of item
            calc_info = analyzer_list[[nr_analyzer]]@results[[nr_results]]$information$calc_info[nr_data], #calc info of item
            plots =lapply(seq_along(analyzer_list[[nr_analyzer]]@results[[nr_results]]$plots),function(nr_combn1) 
              lapply(seq_along(analyzer_list[[nr_analyzer]]@results[[nr_results]]$plots[[nr_combn1]]), function(nr_combn2) 
                lapply(seq_along(analyzer_list[[nr_analyzer]]@results[[nr_results]]$plots[[nr_combn1]]), function(nr_combn2) 
                analyzer_list[[nr_analyzer]]@results[[nr_results]]$plots[[nr_combn1]][[nr_combn2]][[nr_data]]
              )
              ) %>% Reduce(c,.)
            )%>% Reduce(c,.)
          )
        }) %>% 
          setNames(names(analyzer_list[[nr_analyzer]]@results[[nr_results]]$plot_data))
      )
    }) %>% 
      setNames(names(analyzer_list[[nr_analyzer]]@results)) 
  }
  
  # header ----
  # write the header of the rmd file
  write_report_header <- function(title, type="html"){
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
  out_header <- write_report_header(title=ifelse(is.null(title),paste0("metime",Sys.Date()),tools::file_path_sans_ext(file)),
                                    type=device)
  
  # CSS block 
  out_css <- ifelse(
    device == "html",
    "```{css, echo=FALSE}
     #content{max-width:100%}
     #sidebar{bacgkround:black}
     #postamble{border-top:solid 10px black}
     .container{width:100%}
     .button_fun{padding: 0 2px; border:transparent; background:transparent; color:#9F2042; }
     .button_fun:hover{padding: 0 2px; border:transparent; background:transparent; color:#9F2042;font-weight:bold}
     .body_fun{border-radius: 10px; background:rgb(214, 214, 214); padding:10px; margin:10px 0;position:relative}
     .button_copy{margin:10px;padding:5px;border-radius:5px;background:white;color:#9F2042;border-color:transparent;position:absolute;top:0;right:0}
     .button_copy:hover{margin:10px;padding:5px;background:#9F2042;color:white;}
```",
    ""
  )
  
  # JS block
  out_js <- ''
  
  # body -----
  # write the body of the rmd file
  
  ## functions needed to create the body
  ### write information as a text block
  write_report_text<-function(..., type="html"){
    out <- lapply(list(...),function(x)paste0('\n',x,'\n')) %>% 
      unlist()
    return(out)
  }
  
  ### write information as a heading with the option of setting tabs
  write_report_heading <- function(text, level=1, tabset=F, type="html"){
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
  
  ### write table into the rmd file
  write_report_table <- function(table_location, caption=NULL, rownames=F, type="html"){
    table_location <- table_location %>% gsub(pattern='"', replacement = "'")
    if(type=="html"){
      out <- paste0('```{r ,echo=F}\nDT::datatable(',table_location,', caption = "',caption,'", rownames=',rownames, ',' ,
                    'extensions = "Buttons",\n options = list(dom = "Blfrtip",\nbuttons = c("csv", "excel")))\n``` ')
    }
    else if(type=="pdf"){
      out <- paste0('```{r ,echo=F}\nknitr::kable(',table_location,', caption = "',caption,'", rownames=',rownames,'")\n```')
    }
    else{
      out <- ""
    }
    return(out)
  }
  
  ### write plot into the rmd file
  write_report_plot <- function(plot_location, title=NULL, caption=NULL, type="html", interactive=F){
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
  
  ### write information on the analyzer object itself
  write_report_analyzer_info <- function(object_location){
    # function that prints all information of the analyzer object in rmarkdown
    out <- paste0("## Analyzer info\n```{r,echo=F, warning=FALSE, message=FALSE}\n",
                  'results <- ',object_location,'@results\n ',
                  'print(',object_location,')',
                  "\n```\n", collapse="")
    return(out)
  }
  
  ### write function information into rmd file
  write_report_function_info <- function(fun, type="html"){
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
        my_fun_db <- get_db[[grep(pattern=paste0(fun,".Rd"), names(get_db), value=T)]]
        out <- capture.output(tools::Rd2HTML(my_fun_db))[-c(1:7)]
=======
    
    # set parameters -------
    ## set name of the file
    out_file <- ifelse(is.null(file), paste0("metime",Sys.Date()), tools::file_path_sans_ext(file))
    
    # process data -------
    ## save analyzer objects as a list
    if(!is.list(object)){ 
      analyzer_list <- list(object)
      analyzer_list_names <- "object1"
    }else{
      analyzer_list <- object
    }
    
    ## reorder the results list and annotate as results_report
    results_list <- list()
    for(nr_analyzer in seq_along(analyzer_list)){
      results_list[[nr_analyzer]] <- lapply(seq_along(analyzer_list[[nr_analyzer]]@results),function(nr_results){
        list(
          functions_applied = analyzer_list[[nr_analyzer]]@results[[nr_results]]$functions_applied,
          items = lapply(seq_along(analyzer_list[[nr_analyzer]]@results[[nr_results]]$plot_data),function(nr_data){
            list(
              h2=names(analyzer_list[[nr_analyzer]]@results[[nr_results]]$plot_data)[nr_data], # name of item
              plot_data = analyzer_list[[nr_analyzer]]@results[[nr_results]]$plot_data[[nr_data]], #plot data of item
              calc_type = analyzer_list[[nr_analyzer]]@results[[nr_results]]$information$calc_type[nr_data], #calc type of item
              calc_info = analyzer_list[[nr_analyzer]]@results[[nr_results]]$information$calc_info[nr_data], #calc info of item
              plots =lapply(seq_along(analyzer_list[[nr_analyzer]]@results[[nr_results]]$plots),function(nr_combn1) 
                lapply(seq_along(analyzer_list[[nr_analyzer]]@results[[nr_results]]$plots[[nr_combn1]]), function(nr_combn2) 
                  analyzer_list[[nr_analyzer]]@results[[nr_results]]$plots[[nr_combn1]][[nr_combn2]][[nr_data]]
                )%>%Reduce(c,.)
              )
            )
          }) %>% 
            setNames(names(analyzer_list[[nr_analyzer]]@results[[nr_results]]$plot_data))
        )
      }) %>% 
        setNames(names(analyzer_list[[nr_analyzer]]@results)) 
    }
    
    # header ----
    # write the header of the rmd file
    write_report_header <- function(title, type="html"){
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
>>>>>>> Stashed changes
      }
      return(out)
    }
<<<<<<< Updated upstream
    out <- out %>% gsub(pattern='"', replacement = "'", fixed = T)
    return(out)
  }
  
  write_report_pipe_information <- function(list, prefix=NULL,type="html"){
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
      
      pipe_copy <- lapply(seq_along(list), function(x){
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
        paste0(names(list)[x],my_params,' %>% <br>')
      })
      if(length(list)>0){
        #remove errors
        pipe_info[[length(pipe_info)]] <- pipe_info[[length(pipe_info)]] %>% gsub(pattern="%>%", replacement = "", fixed=T)
        pipe_copy[[length(pipe_copy)]] <- pipe_copy[[length(pipe_copy)]] %>% gsub(pattern="%>%", replacement = "", fixed=T)
      }
      
      pipe_modal <- lapply(seq_along(list), function(x){
        paste0("```{r ,echo=F}\n",
               'bsplus::bs_modal(id = "',paste0(prefix,x),'",',
               'title = "',names(list)[x],'",',
               'size = "large",',
               'body = htmltools::HTML("',paste0(write_report_function_info(names(list)[x]),collapse=""),'")',
               ')\n``` \n')
      }) %>% unlist()
      
      out <- c(paste0(
        "\n```{r, echo=F}\n",
        "## Pipe information\n",
        'htmltools::HTML("<div class=body_fun> object %>% <br>',
        "<button class=button_copy onclick='copyToClipboard()'>Copy to Clipboard</button>
=======
    out_header <- write_report_header(title=ifelse(is.null(title),paste0("metime",Sys.Date()),tools::file_path_sans_ext(file)),
                                      type=device)
    
    # CSS block 
    out_css <- ifelse(
      device == "html",
      "```{css, echo=FALSE}
     #content{max-width:100%}
     #sidebar{bacgkround:black}
     #postamble{border-top:solid 10px black}
     .container{width:100%}
     .button_fun{padding: 0 2px; border:transparent; background:transparent; color:#9F2042; }
     .button_fun:hover{padding: 0 2px; border:transparent; background:transparent; color:#9F2042;font-weight:bold}
     .body_fun{border-radius: 10px; background:rgb(214, 214, 214); padding:10px; margin:10px 0;position:relative}
     .button_copy{margin:10px;padding:5px;border-radius:5px;background:white;color:#9F2042;border-color:transparent;position:absolute;top:0;right:0}
     .button_copy:hover{margin:10px;padding:5px;background:#9F2042;color:white;}
```",
      ""
    )
    
    # JS block
    out_js <- ''
    
    # body -----
    # write the body of the rmd file
    
    ## functions needed to create the body
    ### write information as a text block
    write_report_text<-function(..., type="html"){
      out <- lapply(list(...),function(x)paste0('\n',x,'\n')) %>% 
        unlist()
      return(out)
    }
    
    ### write information as a heading with the option of setting tabs
    write_report_heading <- function(text, level=1, tabset=F, type="html"){
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
    
    ### write table into the rmd file
    write_report_table <- function(table_location, caption=NULL, rownames=F, type="html"){
      table_location <- table_location %>% gsub(pattern='"', replacement = "'")
      if(type=="html"){
        out <- paste0('```{r ,echo=F}\nDT::datatable(',table_location,', caption = "',caption,'", rownames=',rownames, ',' ,
                      'extensions = "Buttons",\n options = list(dom = "Blfrtip",\nbuttons = c("csv", "excel")))\n``` ')
      }
      else if(type=="pdf"){
        out <- paste0('```{r ,echo=F}\nknitr::kable(',table_location,', caption = "',caption,'", rownames=',rownames,'")\n```')
      }
      else{
        out <- ""
      }
      return(out)
    }
    
    ### write plot into the rmd file
    write_report_plot <- function(plot_location, title=NULL, caption=NULL, type="html", interactive=F){
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
    
    ### write information on the analyzer object itself
    write_report_analyzer_info <- function(object_location){
      # function that prints all information of the analyzer object in rmarkdown
      out <- paste0("## Analyzer info\n```{r,echo=F, warning=FALSE, message=FALSE}\n",
                    'results <- ',object_location,'@results\n ',
                    'print(',object_location,')',
                    "\n```\n", collapse="")
      return(out)
    }
    
    ### write function information into rmd file
    write_report_function_info <- function(fun, type="html"){
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
          my_fun_db <- get_db[[grep(pattern=paste0(fun,".Rd"), names(get_db), value=T)]]
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
    
    write_report_pipe_information <- function(list, prefix=NULL,type="html"){
      if(type=="html"){
        #gather pipe information
        functions <- list
        names(functions) <- strsplit(functions, split='[(]') %>% lapply(function(x) return(x[1]))
        pipe_info <- lapply(seq_along(functions), function(ind) {
              functions[ind] <- functions[ind] %>% gsub(pattern='"', replacement="'")
              paste0('   <button class=button_fun data-toggle=modal data-target=',paste0("#",prefix,ind),' >', names(functions)[ind],'</button><span>',functions[ind] %>% gsub(pattern=names(functions)[ind], replacement=""),' %>% </span><br>')
          })

        # pipe_info <- lapply(seq_along(list), function(x){
        #   this_info <- data.frame(var = names(list[[x]]) %>% gsub(pattern="[0-9]", replacement = ""), par=as.character(list[[x]]), stringsAsFactors = F)
        #   my_params <-paste0("(",
        #                      lapply(unique(this_info$var), function(y){
        #                        if(y!="mutations"){
        #                          # if mutations then this item was used in mod_mutate or mod_filter
        #                          paste0(y, " = ", ifelse(length(which(this_info$var==y))==1,
        #                                                  paste0("'",this_info$par[which(this_info$var==y)],"'"), 
        #                                                  paste0("c('",paste0(this_info$par[which(this_info$var==y)], collapse="', '"),"')")))
        #                        }else{
        #                          paste0( ifelse(length(which(this_info$var==y))==1,
        #                                         paste0(this_info$par[which(this_info$var==y)],"'"), 
        #                                         paste0("c('",paste0(this_info$par[which(this_info$var==y)], collapse=" , "),")")))
        #                        }
        #                      })%>%
        #                        paste0(collapse = ", "), ")") %>% gsub(pattern='"',replacement = "'",fixed=T)
        #   paste0('   <button class=button_fun data-toggle=modal data-target=',paste0("#",prefix,x),' >', names(list)[x],'</button><span>',my_params,' %>% </span><br>')
        # })
        pipe_info[[length(pipe_info)]] <- pipe_info[[length(pipe_info)]] %>% gsub(pattern="%>%", replacement = "", fixed=T)
        
        # pipe_copy <- lapply(seq_along(list), function(x){
        #   this_info <- data.frame(var = names(list[[x]]) %>% gsub(pattern="[0-9]", replacement = ""), par=as.character(list[[x]]), stringsAsFactors = F)
        #   my_params <-paste0("(",
        #                      lapply(unique(this_info$var), function(y){
        #                        if(y!="mutations"){
        #                          # if mutations then this item was used in mod_mutate or mod_filter
        #                          paste0(y, " = ", ifelse(length(which(this_info$var==y))==1,
        #                                                  paste0("'",this_info$par[which(this_info$var==y)],"'"), 
        #                                                  paste0("c('",paste0(this_info$par[which(this_info$var==y)], collapse="', '"),"')")))
        #                        }else{
        #                          paste0( ifelse(length(which(this_info$var==y))==1,
        #                                         paste0(this_info$par[which(this_info$var==y)],"'"), 
        #                                         paste0("c('",paste0(this_info$par[which(this_info$var==y)], collapse=" , "),")")))
        #                        }
        #                      })%>%
        #                        paste0(collapse = ", "), ")") %>% gsub(pattern='"',replacement = "'",fixed=T)
        #   paste0(names(list)[x],my_params,' %>% <br>')
        # })
        pipe_copy <- lapply(seq_along(functions), function(ind) {
              functions[ind] <- functions[ind] %>% gsub(pattern='"', replacement="'", fixed=TRUE)
              paste0(functions[ind], ' %>% <br>')
          })
        pipe_copy[[length(pipe_copy)]] <- pipe_copy[[length(pipe_copy)]] %>% gsub(pattern="%>%", replacement = "", fixed=T)
        
        pipe_modal <- lapply(seq_along(functions), function(x){
          paste0("```{r ,echo=F}\n",
                 'bsplus::bs_modal(id = "',paste0(prefix,x),'",',
                 'title = "',names(functions)[x],'",',
                 'size = "large",',
                 'body = htmltools::HTML("',paste0(write_report_function_info(names(functions)[x]),collapse=""),'")',
                 ')\n``` \n')
        }) %>% unlist()
        
        out <- c(paste0(
          "\n```{r, echo=F}\n",
          "## Pipe information\n",
          'htmltools::HTML("<div class=body_fun> object %>% <br>',
          "<button class=button_copy onclick='copyToClipboard()'>Copy to Clipboard</button>
>>>>>>> Stashed changes
                      <script>\n function copyToClipboard() {var text = '", paste0(unlist(pipe_copy), collapse=""),"';
var tempInput = document.createElement('input');
            tempInput.value = text;
            document.body.appendChild(tempInput);
            tempInput.select();
            document.execCommand('copy');
            document.body.removeChild(tempInput);
    }\n
  </script>\n",
<<<<<<< Updated upstream
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
      if(length(list)>0){
        out[[length(out)]] <- out[[length(out)]] %>% gsub(pattern="*%>%*", replacement = "", fixed=T)
      }
      out <- paste0("\n > ###### object *%>%* \n", paste0(out, collapse=""))
    }else{
      out = ""
    }
    return(out)
  }
  
  ## write body
  out_body <- lapply(seq_along(results_list), function(nr_analyzer){
    list(
      this_subitem = lapply(seq_along(results_list[[nr_analyzer]]),function(nr_subitem){
        list(
          section_title = write_report_heading(text=names(results_list[[nr_analyzer]])[nr_subitem], tabset=F, level=1, type=device),
          section_analyzer = write_report_analyzer_info(paste0("analyzer_list[[",nr_analyzer,"]]")),
          section_pipe = write_report_pipe_information(results_list[[nr_analyzer]][[nr_subitem]][["functions_applied"]], prefix=nr_subitem, type=device),
          section_plot_info = write_report_heading(text="Plots",level=2, tabset=T,type=device),
          section_body = lapply(seq_along(results_list[[nr_analyzer]][[nr_subitem]]$items),function(nr_plot){
            list(
              plot_heading=write_report_heading(text=names(results_list[[nr_analyzer]][[nr_subitem]]$items)[nr_plot], level=3, tabset=F, type=device),
              plot_info = write_report_text(results_list[[nr_analyzer]][[nr_subitem]]$items[[nr_plot]]$calc_info,type=device),
              plot = lapply(seq_along(results_list[[nr_analyzer]][[nr_subitem]]$items[[nr_plot]]$plots), function(nr_subplot){
                if(is.list(results_list[[nr_analyzer]][[nr_subitem]]$items[[nr_plot]]$plots[[nr_subplot]])){
                  # if it is a list, display as a figure
                  write_report_plot(
                    plot_location = paste0("results_list[[",nr_analyzer,"]][[",nr_subitem,"]]$items[[",nr_plot,"]]$plots[[",nr_subplot,"]]"),
                    type=device, 
                    interactive = interactive)
                  
                }else{
                  # if the item is a data.frame display as table
                  write_report_table(
                    table_location = paste0("results_list[[",nr_analyzer,"]][[",nr_subitem,"]]$items[[",nr_plot,"]]$plots[[",nr_subplot,"]]"),
                    type=device, 
                    interactive = interactive)
                }
              }) %>% unlist()
            )
          }) %>% unlist()
        )
      }) %>% unlist()
    )
  }) %>% unlist()
  
  # Footer ####
  out_footer <- "# Session info\n```{r session_info,include=T, warning=FALSE, message=FALSE}\nprint(sessionInfo())\n```\n"
  
  # Write xlsx
  if(write_results){
    write_results(object = object, file=out_file)
  }
  
  # write report ####
  writeLines(
    text=c(
      out_header,
      out_css,
      out_js, 
      out_body,
      out_footer), # html code
    paste0(out_file,".rmd") #path for saving the file
  )
  
  # convert rmd to html
  output_format = NULL
  if(device %in% c("pdf","word")) output_format <- paste0(device, "_document")
  rmarkdown::render(input=paste0(out_file,".rmd"),output_format = output_format)
})
    
=======
          '<div style= margin-left:5px>',
          paste0(pipe_info, collapse=""), 
          '</div>',
          '</div>")\n```\n'
        ),
        pipe_modal
        )
      }
      else if(type=="pdf"){
        out <- lapply(seq_along(functions), function(x){
          # this_info <- data.frame(var = names(list[[x]]) %>% gsub(pattern="[0-9]", replacement = ""), par=as.character(list[[x]]), stringsAsFactors = F)
          # my_params <-paste0("(",lapply(unique(this_info$var), function(y) 
          #   paste0(y, " = ", ifelse(length(which(this_info$var==y))==1,
          #                           paste0("'",this_info$par[which(this_info$var==y)],"'"), 
          #                           paste0("c('",paste0(this_info$par[which(this_info$var==y)], collapse="', '"),"')")))) %>%
          #     paste0(collapse = ", "), ")") %>% gsub(pattern='"',replacement = "'",fixed=T)
          out <- paste0("\n >> ###### ", functions[x], " *%>%*\n")
        })
        out[[length(out)]] <- out[[length(out)]] %>% gsub(pattern="*%>%*", replacement = "", fixed=T)
        out <- paste0("\n > ###### object *%>%* \n", paste0(out, collapse=""))
      }else{
        out = ""
      }
      return(out)
    }
    
    ## write body
    out_body <- lapply(seq_along(results_list), function(nr_analyzer){
      list(
        this_subitem = lapply(seq_along(results_list[[nr_analyzer]]),function(nr_subitem){
          list(
            section_title = write_report_heading(text=names(results_list[[nr_analyzer]])[nr_subitem], tabset=F, level=1, type=device),
            section_analyzer = write_report_analyzer_info(paste0("analyzer_list[[",nr_analyzer,"]]")),
            section_pipe = write_report_pipe_information(results_list[[nr_analyzer]][[nr_subitem]][["functions_applied"]], prefix=nr_subitem, type=device),
            section_plot_info = write_report_heading(text="Plot",level=2, tabset=T,type=device),
            section_body = lapply(seq_along(results_list[[nr_analyzer]][[nr_subitem]]$items),function(nr_plot){
              list(
                plot_heading=write_report_heading(text=names(results_list[[nr_analyzer]][[nr_subitem]]$items)[nr_plot], level=3, tabset=F, type=device),
                plot_info = write_report_text(results_list[[nr_analyzer]][[nr_subitem]]$items[[nr_plot]]$calc_info,type=device),
                plot = lapply(seq_along(results_list[[nr_analyzer]][[nr_subitem]]$items[[nr_plot]]$plots), function(nr_subplot){
                  if(is.list(results_list[[nr_analyzer]][[nr_subitem]]$items[[nr_plot]]$plots[[nr_subplot]])){
                    # if it is a list, display as a figure
                    write_report_plot(
                      plot_location = paste0("results_list[[",nr_analyzer,"]][[",nr_subitem,"]]$items[[",nr_plot,"]]$plots[[",nr_subplot,"]]"),
                      type=device, 
                      interactive = interactive)
                    
                  }else{
                    # if the item is a data.frame display as table
                    write_report_table(
                      table_location = paste0("results_list[[",nr_analyzer,"]][[",nr_subitem,"]]$items[[",nr_plot,"]]$plots[[",nr_subplot,"]]"),
                      type=device, 
                      interactive = interactive)
                  }
                }) %>% unlist()
              )
            }) %>% unlist()
          )
        }) %>% unlist()
      )
    }) %>% unlist()
    
    # Footer ####
    out_footer <- "# Session info\n```{r session_info,include=T, warning=FALSE, message=FALSE}\nprint(sessionInfo())\n```\n"
    
    # Write xlsx
    if(write_results){
      write_results(object = object, results_index=NULL)
    }
    # write report ####
    writeLines(
      text=c(
        out_header,
        out_css,
        out_js, 
        out_body,
        out_footer), # html code
      paste0(out_file,".rmd") #path for saving the file
    )
    
    # convert rmd to html
    output_format = NULL
    if(device %in% c("pdf","word")) output_format <- paste0(device, "_document")
    rmarkdown::render(input=paste0(out_file,".rmd"),output_format = output_format)
  })
>>>>>>> Stashed changes
