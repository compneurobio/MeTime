#' Check normality of metabolites and append results to col_data
#' @description A method applied to data within the S4 object of class "metime_analyser" to check normality of the metabolites
#' and generate a results section for easy visualization,
#' @param object a S4 object of the class "metime_analyzer".
#' @param which_data a character to define which dataset is to be used.
#' @param method a character vector specifying the method to check for normality. By default all available options are chosen: c("shapiro", "ks")". See details for more information.
#' @param cols_for_meta a list of named character vectors for obtaining metainformation of metabolites
#' @param name name of the results. Default is set to "calc_col_stats_1"
#' @details Shapiro: Performs the Shapiro-Wilk test of normality (normal distribution if shapiro_pval > 0.05). This test is based on correlation between the data and the corresponding normal scores.
#' Ks: Performs a one-sample Kolmogorov-Smirnov test (normal distribution if ks_pval > 0.05). The ks.test compares the metabolite values against pnorm.
#' @seealso [shapiro.test] [ks.test] [get_metadata_for_columns]
#' @return a  S4 object of class "metime_analyser" with results appended to the col_data of which_data
#' @export
setGeneric("calc_col_stats", function(object, which_data, method=c("shapiro", "ks"), cols_for_meta, name="calc_col_stats_1") standardGeneric("calc_col_stats"))
setMethod("calc_col_stats", "metime_analyser", function(object, which_data, method=c("shapiro", "ks"), cols_for_meta, name="calc_col_stats_1") {
  out <- object
  # sanity checks
  if(!all(which_data %in% names(object@list_of_data))) {
      warning("add_col_stats(): which_data not in metime_analyzer, or more than one which_data selected")
  } else if(!all(method %in% c("shapiro", "ks"))) {
    warning("add_col_stats(): type can only be 'shapiro', 'ks'")
  } else {
    # calculate iterations across which_data
    if(grep(name, names(object@results)) %>% length() >=1) {
                        warning("name of the results was previously used, using a different name")
                        index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
                        index <- c(0:9)[grep(index, 0:9)+1]
                        name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
                }
    data_list <- lapply(which_data, function(i) {
      # iterate per which_data column
      lapply(out@list_of_col_data[[i]]$id, function(x) {
        this_out <- data.frame(id=x, stringsAsFactors = F)
        if("shapiro" %in% method){
          # shapiro test from base::sharpio.test
          my_shapiro <- shapiro.test(object@list_of_data[[i]][, x])
          this_out <- this_out %>% 
            dplyr::mutate(id=x,
                          shapiro_pval=as.numeric(my_shapiro$p.value),
                          shapiro_statistic=as.numeric(my_shapiro$statistic),
                          shapiro_normal = ifelse(as.numeric(my_shapiro$p.value)>0.05, TRUE,FALSE))
        }
        if("ks" %in% method){
          # ks test from base::ks.test
          my_ks <- ks.test(x=object@list_of_data[[i]][, x], y='pnorm')
          this_out <- this_out %>% 
            dplyr::mutate(id=x,
                          ks_pval=as.numeric(my_ks$p.value),
                          ks_statistic=as.numeric(my_ks$statistic),
                          ks_normal = ifelse(as.numeric(my_ks$p.value)>0.05, TRUE,FALSE))
        }
        this_out
      }) %>% 
        plyr::rbind.fill()

    }) 
    data_list <- setNames(data_list, which_data)
    data_list <- lapply(data_list, function(df) {
      rownames(df) <- df$id
      df
    })
    if(is.null(cols_for_meta)) {
      metadata <- NULL
    } else {
      metadata <- get_metadata_for_columns(object=object, 
                                which_data=which_data, columns=cols_for_meta)
    }
    col_stats <- if (length(which_data) == 1) data_list[[1]] else data_list
    out <- get_make_results(object=object, data=col_stats, metadata=metadata, calc_type="normality", 
                  calc_info=paste("Normality statistics of ", which_data, sep=""), name=name) %>%
           add_function_info(function_name="calc_col_stats", 
                  params=list(which_data=which_data, method=method, cols_for_meta=cols_for_meta))
  }
    

    return(out)
    
  })
