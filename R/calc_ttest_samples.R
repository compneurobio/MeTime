
#' Function to calculate students t-test between samples at different timepoints
#' @description Method for S4 object of class metime_analyser for performing t-test
#' @param object S4 object of class metime_analyser
#' @param which_data dataset or datasets to be used for the analysis
#' @param timepoints timepoints of interest to perform the test on
#' @param type type of ttest to be used either "two.sided", "less", or "greater"
#' @param paired Logical to perform paired t.test or not
#' @return plotter object with t-test results
#' @export
setGeneric("calc_ttest_samples", function(object, which_data, timepoints, type, paired=TRUE) standardGeneric("calc_ttest_samples"))
setMethod("calc_ttest_samples", "metime_analyser", function(object, which_data, timepoints, type, paired=TRUE) {
        if(length(which_data) > 1) {
          object <- mod_extract_common_samples(object)
          combinations <- as.data.frame(t(combn(timepoints, 2)))
          list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
          data <- do.call(unname(list_of_data), cbind)
          timepoints <- as.character(unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[2]))))
          samples <- as.character(unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[1]))))
          data <- as.data.frame(cbind(data, timepoints, samples))
        } else {
          data <- as.data.frame(object@list_of_data[[which_data]])
          combinations <- as.data.frame(t(combn(timepoints, 2)))
          timepoints <- as.character(unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[2]))))
          samples <- as.character(unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[1]))))
          data <- as.data.frame(cbind(data, timepoints, samples))
          name_combs <- cbind(unique(data$samples), unique(data$samples))
          colnames(name_combs) <- c("s1", "s2")
        }
        out <- lapply(1:nrow(combinations), function(x) {
              t1 <- as.character(combinations[x, 1])
              t2 <- as.character(combinations[x, 2])
              t1_data <- data %>% filter(timepoints==t1) 
              t2_data <- data %>% filter(timepoints==t2)
              rownames(t1_data) <- t1_data$samples
              rownames(t2_data) <- t2_data$samples
              t1_data <- t1_data %>% dplyr::select(-samples, -timepoints) %>% t()
              t2_data <- t2_data %>% dplyr::select(-samples, -timepoints) %>% t()
              common_samples <- intersect(colnames(t1_data), colnames(t2_data))
              name_combs <- cbind(common_samples, common_samples)
              out <- parallel::mclapply(name_combs, function(x) {
                  s1 <- t1_data[ ,x$V1]
                  s2 <- t2_data[ ,x$V2]
                  result <- rstatix::t_test(s1 ~ s2, paired=paired, alternative=type)
                  return(result)
                }, max.cores=4) %>% do.call(what=rbind.data.frame)
              return(results)
          })
        return(out)
  })


