
#Metabolites with pvalues, fold change(statistic) and timepoints information

#' An automated fucntion to calculate GGM from genenet longitudnal version
#' @description automated funtion that can be applied on metime_analyser object to obtain geneNet network along with threshold used
#' @param object S4 object of cĺass metime_analyser
#' @param which_data a character or a character vector naming the datasets of interest
#' @param threshold type of threshold to be used for extracting significant edge
#' @param name character vector for naming the results
#' @param stratifications List to stratify data into a subset. Usage list(name=value)
#' @param all Logical to get all edges without any cutoff.
#' @param cols_for_meta a list of character vectors for getting metadata for columns for plotting purposes
#' @param covariates covariates to be used for this analysis
#' @param ... addtional arguments for genenet network
#' @return Network data as a plotter object 
#' @export
setGeneric("calc_ggm_genenet_longitudnal", function(object, which_data, threshold, all, cols_for_meta, covariates, stratifications, name, ...) standardGeneric("calc_ggm_genenet_longitudnal"))
setMethod("calc_ggm_genenet_longitudnal", "metime_analyser", function(object, which_data, threshold, all, cols_for_meta, covariates, stratifications, name, ...) {
    #sanity checks
    stopifnot(time %in% names(stratifications))
    stopifnot(threshold %in% c("li","bonferroni","FDR", NULL))
    #Extracting data that is needed
    object@list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
    object@list_of_row_data <- object@list_of_row_data[names(object@list_of_row_data) %in% which_data]
    #converting prepped data to full ggm network
    if(length(which_data > 1)) {
        object <- mod_extract_common_samples(object)
        list_of_data <- unname(object@list_of_data)
        data <- do.call(cbind, list_of_data)
    } else {
        data <- object@list_of_data[[1]]
    } 
    data <- na.omit(data)
    data$subject <- rownames(data) %>% gsub(pattern="_[a-z|A-Z][0-9]+", replacement="")
    data$time <- rownames(data) %>% gsub(pattern="[a-z|A-Z][0-9]+_", replacement="")
    my_subjects <- data %>%     
              dplyr::select(all_of(c("subject","time"))) %>% 
              dplyr::filter(get("time") %in% stratifications$time)%>% 
              dplyr::group_by(subject) %>% 
              dplyr::count(subject) %>% 
              dplyr::filter(n == length(stratifications$time))
                  
    data <- data %>% 
            dplyr::filter(time %in% timepoints,
            subject %in% my_subjects[["subject"]])
    rm_col = intersect(names(data), c("adni_id","RID","rid","time","tp","subject", "id"))
    vars <- data %>% select(-c(all_of(rm_col))) %>% names()


    # get full data
    data <- data %>% dplyr::arrange(time, subject) 
    n_subject = unique(data$subject) %>% length() %>% as.numeric()
    name_tp = as.numeric(unlist(lapply(strsplit(unique(data$time), split="t"), function(x) return(x[2]))))

    data <- longitudinal::as.longitudinal(x=as.matrix(data[,vars]), repeats=n_subject, time=name_tp)

    network <- get_ggm_genenet(data=data, threshold=threshold, all=all, ...)
    network <- network[!network$node1 %in% covariates, ]
    network <- network[!network$node2 %in% covariates, ]
    metadata <- get_metadata_for_columns(object=object, which_data=which_data, columns=cols_for_meta, 
                 names=c("name", "pathway"), index_of_names="id")
    out <- get_results(object=object, data=network, metadata=metadata, calc_type="genenet_ggm", 
      calc_info=paste("Longitudinal GeneNet GGM results for:", paste(which_data, collapse=" & "), sep=" "),
      name=name)
       
    return(out)
})


