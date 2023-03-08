
#' Function to calculate metabotype conservation index
#' @description Method applied on the object metime_analyser to calculate the metabotype conservation index
#' @examples #calculating metabotype_conservation_index 
#' out <- calc_metabotype_conservation(object=metime_analyser_object, which_data="Name of the dataset")
#' @param object An object of class metime_analyser
#' @param which_data Name of the dataset to be used
#' @param verbose Information provided on steps being processed
#' @param cols_for_meta Character vector to define column names that are to be used for plotting purposes
#' @param name character vector to define the results. Should be equal to length of which_data
#' @param stratifications List to stratify data into a subset. Usage list(name=value)
#' @importClassesFrom metime_analyser
#' @return List of conservation index results
#' @export
setGeneric("calc_conservation_metabotype", function(object, which_data, verbose=F, cols_for_meta, stratifications, name="calc_conservation_metabotype_1") standardGeneric("calc_conservation_metabotype"))
setMethod("calc_conservation_metabotype", "metime_analyser", function(object, which_data, verbose=F, cols_for_meta, stratifications, name="calc_conservation_metabotype_1") {

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
    data_list <- get_stratified_data(which_data=which_data, object=object, 
      stratifications=stratifications)
    # extract data
    data <- data_list[["data"]]
    row_data <- data_list[["row_data"]]
    this_data <- data %>% 
      dplyr::mutate(id = rownames(.[])) %>% 
      dplyr::left_join(row_data[,index_var], by = "id") %>% 
      dplyr::mutate(time = gsub(x=id, pattern="[a-z|A-Z][-|0-9]+_", replacement="")) %>% 
      `rownames<-`(.[,"id"]) %>% 
      dplyr::arrange(subject, time)
    # calculate correlations
    data_cor <- this_data %>%
      dplyr::select(-id, -time, -subject) %>% 
      t() %>% 
      data.frame() %>% 
      cor(use = "pairwise.complete.obs")
    
    # index time point combinations
    index_time_combinations <- utils::combn(x = unique(this_data$time), 
                                            m = 2, simplify = T)
    # get all TP combinations in a prelim. results
    out_sum <- lapply(1:ncol(index_time_combinations), function(x) {
       this_out <- dplyr::full_join(x=this_data[this_data$time == index_time_combinations[1,x],c("id","subject","time")],
                                   y=this_data[this_data$time == index_time_combinations[2,x],c("id","subject","time")],
                                   by="subject") %>% 
         dplyr::rename("id_from"="id.x","id_to"="id.y", "time_from"="time.x","time_to"="time.y") %>%
          na.omit() # Fixed bug below by adding this line
       #add results 
       this_out$cor <- data_cor[this_out$id_from,this_out$id_to] %>% diag() # Throws error
       this_out$rank <- as.data.frame(data_cor[this_out$id_from,this_out$id_to]* -1)%>% 
         dplyr::mutate_all(~ rank(.x,ties.method = "min")) %>%
         as.matrix() %>% 
         diag()
       
       this_out %>% 
         dplyr::mutate(ci = (1 - ((rank - 1)/(nrow(.[]) - 1)))) %>% 
         dplyr::arrange(desc(ci)) %>% 
         dplyr::mutate(subject_rank = nrow(.[]):1) %>%
         dplyr::mutate(x=subject_rank,
                       id=subject,
                       nsubject=nrow(.[])) %>% 
         dplyr::select(x, ci, id, time_from, time_to, nsubject, rank, cor, id_from, id_to) %>% 
         `rownames<-`(.[,"id_from"]) 
    })
    
    metadata <- get_metadata_for_rows(object = object, 
                                      which_data = i, 
                                      columns = cols_for_meta)
    out <- list()
    combinations <- lapply(1:ncol(index_time_combinations), function(y) {
              t <- paste(index_time_combinations[,y], collapse="vs")
              return(t)
      }) %>% unlist()

    out <- get_make_results(object=object, data = out_sum, 
                              metadata = metadata, 
                              calc_type = rep("CI_metabotype", each=length(out_sum)), 
                              calc_info = paste("metabotype_CI_", i, "_", combinations, sep = ""),
                              name=name)
    out <- add_function_info(object=out, function_name="calc_conservation_metabotype", 
        params=list(which_data=which_data, verbose=verbose, cols_for_meta=cols_for_meta, 
            name=name[i], stratifications=stratifications))
  }
  return(out)
})

#object <- calc_conservation_metabotype(object=object, which_data="lipid_data", 
#stratifications=list(time=c("t0", "t12", "t24")), verbose=TRUE, 
#cols_for_meta=c("ADNI_MEM", "ADNI_LAN", "ADNI_EF", "APOEGrp",  "DXGrp_longi", "PTGENDER", "Age", "BMI"), name="CI_test")

