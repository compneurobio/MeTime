

#' Function to calculate colinearity 
#' @description Function to calculate colinearity in a dataset
#' @param object An S4 object of class metime_analyser
#' @param which_data Dataset to check for colinearity
#' @param cols_for_meta character vector to columns needed in metadata. Id and class name is needed
#' @param show_all logical. True will only filter out colinear data
#' @return plotter object with data for heatmap information
#' @export
setGeneric("calc_colinearity", function(object, which_data, cols_for_meta, show_all) standardGeneric("calc_colinearity"))
setMethod("calc_colinearity", "metime_analyser", function(object, which_data, cols_for_meta, show_all) {
      stopifnot(which_data %in% names(object@list_of_data))
      stopifnot(length(names(object@list_of_data[[which_data]]))==length(unique(names(object@list_of_data[[which_data]]))))
      my_combn <- combn(names(object@list_of_data[[which_data]]),2) %>% t() %>% as.data.frame() # setup combinations for later pairwise calculation

      out <- lapply(1:nrow(my_combn), function(i) {
          my_x=as.numeric(object@list_of_data[[which_data]][,my_combn$V1[i]])
          my_y=as.numeric(object@list_of_data[[which_data]][,my_combn$V2[i]])
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
      metadata <- get_metadata_for_columns(object=object, which_data=which_data, columns=list(cols_for_meta), 
                      names=c("id", "class"), 
                      index_of_names="id")
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
          out <- get_make_plotter_object(data=out, metadata=NULL, calc_type="colinearity", 
                  calc_info=paste("colinearity test of: ", which_data, sep=""),
                  plot_type="heatmap", style="ggplot")
          return(out)
      } else {
          out <- out[out$colinear==TRUE, ]
          if(dim(out)[1]==0) {
              cat("There are no colinear vectors in the data")
              cat("\n")
              cat("Exiting without making a plotter object")
              return(out)
          } else {
              for(i in 1:length(out$med_class_1)) {
                  out$med_class_1[i] <- as.character(metadata[metadata$id %in% out$med_class_1[i], "class"])
                  out$med_class_2[i] <- as.character(metadata[metadata$id %in% out$med_class_2[i], "class"])
              }
              out <- get_make_plotter_object(data=out, metadata=NULL, calc_type="colinearity", 
                  calc_info=paste("colinearity test of: ", which_data, sep=""),
                  plot_type="heatmap", style="ggplot")
              return(out)
          }
      }
  })


