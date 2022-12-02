
#' Function to calculate metabolite conservation index
#' @description Method applied on the object metime_analyser to calculate the metabotype conservation index
#' @examples #calculating metabolite_conservation_index 
#' out <- calc_metabolite_conservation(object=metime_analyser_object, which_data="Name of the dataset")
#' @param object An object of class metime_analyser
#' @param which_data Name of the dataset to be used
#' @param timepoints character vector with timepoints of interest
#' @param verbose Information provided on steps being processed
#' @param cols_for_meta A list of a Character vector to define column names that are to be used for plotting purposes
#' @return List of conservation index results
#' @export
setGeneric("calc_conservation_metabolite", function(object, which_data, timepoints, verbose, cols_for_meta) standardGeneric("calc_conservation_metabolite"))
setMethod("calc_conservation_metabolite", "metime_analyser", function(object, which_data, timepoints, verbose=F, cols_for_meta) {
  #define data to be processed
  #data_position <- which(names(object@list_of_data) %in% which_data)
  out=list()
  for (i in which_data) {
    
    # index variables to use
    index_var <- c("id","time","subject")
    
    # extract and reshape data
    this_data <- object@list_of_data[[which(names(object@list_of_data)==i)]] %>% 
      dplyr::mutate(id = rownames(.[])) %>% 
      dplyr::left_join(object@list_of_row_data[[which(names(object@list_of_data)==i)]][,index_var], by = "id") %>% 
      dplyr::mutate(time = gsub(x=time, pattern="t", replacement="")) %>% 
      `rownames<-`(.[,"id"]) %>% 
      dplyr::arrange(subject, time)
  
    this_data_reshaped <- lapply(unique(this_data$time),function(x){
      this_data %>% 
        dplyr::filter(time==x) %>%
        dplyr::select(-time,-id) %>% 
        `colnames<-`(ifelse(colnames(.[])=="subject","subject", paste0(colnames(.[]) ,"timepoint",x)))
    }) %>% 
      plyr::join_all(by="subject")
    
    # calculate correlations
    data_cor <- this_data_reshaped %>%
      dplyr::select(-subject) %>% 
      cor(use = "pairwise.complete.obs")
    
    # index time point combinations
    index_time_combinations <- utils::combn(x = unique(this_data$time), 
                                            m = 2, simplify = T)
    # get all TP combinations in a prelim. results
    out_sum <- lapply(1:ncol(index_time_combinations), function(x){
      this_out <- data.frame(
        id = setdiff(names(this_data), index_var),
        time_from=index_time_combinations[1,x],
        time_to=index_time_combinations[2,x],
        stringsAsFactors = F
      ) %>% 
        dplyr::mutate(id_from=paste0(id, "timepoint",time_from),
                      id_to=paste0(id, "timepoint",time_to)
                      )
      #add results 
      this_out$cor <- data_cor[this_out$id_from,this_out$id_to] %>% diag()
      this_out$rank <- as.data.frame(data_cor[this_out$id_from,this_out$id_to]* -1)%>% 
        dplyr::mutate_all(~ rank(.x,ties.method = "min")) %>%
        as.matrix() %>% 
        diag()
      
      this_out %>% 
        dplyr::mutate(ci = (1 - ((rank - 1)/(nrow(.[]) - 1)))) %>% 
        dplyr::arrange(desc(ci)) %>% 
        dplyr::mutate(x=nrow(.[]):1,
                      y=ci, 
                      n=nrow(.[])) %>% 
        dplyr::select(x,y, ci, id, time_from, time_to, n, rank, cor, id_from,id_to) %>% 
        `rownames<-`(.[,"id_from"])
    })
    
    metadata <- get_metadata_for_columns(object = object, 
                                         which_data = i, 
                                         columns = cols_for_meta, 
                                         names = c("name", "group"), 
                                         index_of_names = "id")
    out <- list()
    out[[i]]<- lapply(1:ncol(index_time_combinations), function(x) {
      get_make_plotter_object(data = out_sum[[x]], 
                              metadata = metadata, 
                              calc_type = "CI", 
                              calc_info = paste("metabotype_CI_", i, "_", x, sep = ""), 
                              plot_type = "dot", 
                              style = "ggplot")
    })
  }
  return(out)
})


