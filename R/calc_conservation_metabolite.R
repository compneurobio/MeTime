
#' Function to calculate metabolite conservation index
#' @description Method applied on the object metime_analyser to calculate the metabotype conservation index
#' @examples #calculating metabolite_conservation_index 
#' out <- calc_metabolite_conservation(object=metime_analyser_object, which_data="Name of the dataset")
#' @param object An object of class metime_analyser
#' @param which_data Name of the dataset to be used
#' @param verbose Information provided on steps being processed
#' @param cols_for_meta A list of a Character vector to define column names that are to be used for plotting purposes
#' the characters should be named the same way eg: list(lipid_data=c(id="id", sub_pathway="sub_pathway"), nmr_data=c(id="id", sub_pathway="Group"))
#' @param name character vector to define the name of the results generated. length should be equal to which_data
#' @param stratifications list to stratify the data used 
#' @return conservation index results that are added to the object
#' @export
setGeneric("calc_conservation_metabolite", function(object, which_data, verbose=F, cols_for_meta=NULL, stratifications=NULL, name="conservation_index_metabolite_1") standardGeneric("calc_conservation_metabolite"))
setMethod("calc_conservation_metabolite", "metime_analyser", function(object, which_data, verbose=F, cols_for_meta=NULL, stratifications=NULL, name="conservation_index_metabolite_1") {
  #define data to be processed
  #data_position <- which(names(object@list_of_data) %in% which_data)
  #test this now
  stopifnot(length(which_data)==length(name))
  if(length(which_data)==1) {
    if(grep(name, names(object@results)) %>% length() >=1) {
      warning("name of the results was previously used, using a different name")
      index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
      index <- c(0:9)[grep(index, 0:9)+1]
      name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
    }
  }
  
  out=list()
  for (i in which_data) {
    
    # index variables to use
    index_var <- c("id","time","subject")
    data_list <- get_stratified_data(which_data=which_data, object=object, stratifications=stratifications)
    data <- data_list[["data"]]
    row_data <- data_list[["row_data"]]
    # extract and reshape data
    this_data <- data %>% 
      dplyr::mutate(id = rownames(.[])) %>% 
      dplyr::left_join(row_data[,index_var], by = "id") %>% 
      dplyr::mutate(time = gsub(x=id, pattern="[a-z|A-Z][-|0-9]+_", replacement="")) %>% 
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
      
      this_out <- this_out %>% 
        dplyr::mutate(ci = (1 - ((rank - 1)/(nrow(.[]) - 1)))) %>% 
        dplyr::arrange(desc(ci)) %>% 
        dplyr::mutate(x=nrow(.[]):1,
                      y=ci, 
                      n=nrow(.[])) %>% 
        dplyr::select(x,y, ci, id, time_from, time_to, n, rank, cor, id_from, id_to) %>% 
        `rownames<-`(.[,"id"])
      return(this_out)
    })

    metadata <- get_metadata_for_columns(object = object, 
                                         which_data = i, 
                                         columns = cols_for_meta, 
                                         names = names(cols_for_meta[[1]]), 
                                         index_of_names = "id")
    out <- list()
    combinations <- lapply(1:ncol(index_time_combinations), function(y) {
              t <- paste(index_time_combinations[,y], collapse="vs")
              return(t)
      }) %>% unlist()

      out <- get_make_results(object=object, data = out_sum, 
                                metadata = metadata, 
                                calc_type = rep("CI_metabolite", each=length(out_sum)), 
                                calc_info = paste("metabolite_CI_", i, "_", combinations, sep = ""),
                                name=name)
      out <- add_function_info(object=out, function_name="calc_conservation_metabolite", 
          params=list(which_data=which_data, verbose=verbose, cols_for_meta=cols_for_meta, 
              name=name[i], stratifications=stratifications))
  }
  return(out)
})


