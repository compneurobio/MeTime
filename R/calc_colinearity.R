

#' Function to calculate colinearity 
#' @description Function to calculate colinearity in a dataset
#' @param object An S4 object of class metime_analyser
#' @param which_data Dataset to check for colinearity
#' @param cols_for_meta list of character vectors of column names needed in metadata. Id and class name is needed
#' @param show_all logical. True will only filter out colinear data
#' @param name character to define the name of the result
#' @param stratifications List to stratify data into a subset. Usage list(name=values)
#' @return plotter object with data for heatmap information
#' @export
setGeneric("calc_colinearity", function(object, which_data, cols_for_meta, show_all=F, name="calc_colinearity_1", stratifications=list()) standardGeneric("calc_colinearity"))
setMethod("calc_colinearity", "metime_analyser", function(object, which_data, cols_for_meta, show_all=F, name="calc_colinearity_1", stratifications=list()) {
      stopifnot(which_data %in% names(object@list_of_data))
      stopifnot(length(names(object@list_of_data[[which_data]]))==length(unique(names(object@list_of_data[[which_data]]))))
      if(grep(name, names(object@results)) %>% length() >=1) {
          warning("name of the results was previously used, using a different name")
          index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
          index <- c(0:9)[grep(index, 0:9)+1]
          name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
      }
      my_combn <- combn(names(object@list_of_data[[which_data]]),2) %>% t() %>% as.data.frame() # setup combinations for later pairwise calculation
      if(length(stratifications)>=1) {
        data <- object@list_of_data[[which_data]]
        row_data <- object@list_of_row_data[[which_data]]
        stratifications <- lapply(names(stratifications), function(x) {
              row_data <- row_data[row_data[,x] %in% stratifications[[x]], ]
              return(stratifications[[x]]) 
          })
        data <- data[rownames(data) %in% rownames(row_data), ]
      } else {
        data <- object@list_of_data[[which_data]]
        row_data <- object@list_of_row_data[[which_data]]
      }
      out <- lapply(1:nrow(my_combn), function(i) {
          my_x=as.numeric(data[,my_combn$V1[i]])
          my_y=as.numeric(data[,my_combn$V2[i]])
          chi_test = stats::chisq.test(x = my_x, 
                                 y = my_y, 
                                 correct=FALSE)
          out = data.frame(
            med_class_1=my_combn$V1[i],
            med_class_2= my_combn$V2[i],
            chi_statistic=chi_test$statistic,
            Cramers_V= chi_test$statistic / (length(my_x) * (min(length(unique(my_x)),length(unique(my_y))) - 1)),
            stringsAsFactors = F
          )
      }) %>%
      do.call(what=rbind.data.frame) %>% 
      dplyr::mutate(colinear=ifelse(Cramers_V>.5, T,F))
      if(is.null(cols_for_meta)) {
        metadata <- NULL
      } else {
        metadata <- get_metadata_for_columns(object=object, 
            which_data=which_data, columns=cols_for_meta)
        metadata <- list(metadata)
      }
      if(show_all) {
          class_1_name <- c()
          class_2_name <- c()
          classes <- c()
          for(i in 1:length(out$med_class_1)) {
              class_1_name[i] <- metadata[metadata$id %in% out$med_class_1[i], "class"]
              class_2_name[i] <- metadata[metadata$id %in% out$med_class_2[i], "class"]
              classes[i] <- paste(as.character(class_1_name[i]), as.character(class_2_name[i]), sep=" - ")
          }
          out <- cbind(out, classes)
          out <- get_make_results(object=object, data=list(out), metadata=metadata, calc_type="colinearity", 
                  calc_info=paste("colinearity test of: ", which_data, sep=""), name=name) %>%
                add_function_info(function_name="calc_colinearity", 
                  params=list(which_data=which_data, cols_for_meta=cols_for_meta, show_all=show_all, stratifications=stratifications))
          return(out)
      } else {
          out <- out[out$colinear==TRUE, ]
          if(dim(out)[1]==0) {
              cat("There are no colinear vectors in the data")
              cat("\n")
              cat("Exiting without updating the object")
              return(object)
          } else {
              for(i in 1:length(out$med_class_1)) {
                  out$med_class_1[i] <- as.character(metadata[metadata$id %in% out$med_class_1[i], "class"])
                  out$med_class_2[i] <- as.character(metadata[metadata$id %in% out$med_class_2[i], "class"])
              }
              out <- get_make_results(object=object, data=list(colinearity=out), metadata=metadata, calc_type="colinearity", 
                  calc_info=paste("colinearity test of: ", which_data, sep=""), name=name)
              out <- out %>% add_function_info(function_name="calc_colinearity", 
                  params=list(which_data=which_data, cols_for_meta=cols_for_meta, show_all=show_all, stratifications=stratifications)) %>%
                  update_plots(type="colinearity")
              return(out)
          }
      }
  })


