

#' Function to calculate students t-test between metabolites at different timepoints
#' @description Method for S4 object of class metime_analyser for performing t-test
#' @param object S4 object of class metime_analyser
#' @param which_data dataset or datasets to be used for the analysis
#' @param timepoints timepoints of interest to perform the test on
#' @param split_var split variable for testing such as diagnostic group etc
#' @param type type of ttest to be used either "two.sided", "less", or "greater"
#' @param paired Logical to perform paired t.test or not
#' @return plotter object with t-test results
#' @export
setGeneric("calc_ttest_metabolites", function(object, which_data, timepoints, split_var, type, paired) standardGeneric("calc_ttest_metabolites"))
setMethod("calc_ttest_metabolites", "metime_analyser", function(object, which_data, timepoints, split_var, type, paired) {
      if(length(which_data) > 1) {
        object <- mod_extract_common_samples(object)
        combinations <- as.data.frame(t(combn(timepoints, 2)))
        list_of_data <- object@list_of_data[names(object@list_of_data) %in% which_data]
        data <- do.call(unname(list_of_data), cbind)
        name_combs <- as.data.frame(cbind(names(data), names(data)))
        colnames(name_combs) <- c("m1", "m2")
        timepoints <- as.character(unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[2]))))
        samples <- as.character(unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[1]))))
        data <- as.data.frame(cbind(data, timepoints, samples))
      } else {
        combinations <- as.data.frame(t(combn(timepoints, 2)))
        data <- object@list_of_data[[which_data]]
        name_combs <- as.data.frame(cbind(names(data), names(data)))
        colnames(name_combs) <- c("m1", "m2")
        timepoints <- as.character(unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[2]))))
        samples <- as.character(unlist(lapply(strsplit(rownames(data), split="_"), function(x) return(x[1]))))
        data <- as.data.frame(cbind(data, timepoints, samples))
      }
      out <- lapply(1:nrow(combinations), function(x) {
            t1 <- as.character(combinations[x, 1])
            t2 <- as.character(combinations[x, 2])
            t1_data <- data %>% filter(timepoints==t1)
            t2_data <- data %>% filter(timepoints==t2)
            if(paired) {
                common_samples <- intersect(t1_data$samples, t2_data$samples)
                t1_data <- t1_data[t1_data$samples %in% common_samples, ]
                t2_data <- t2_data[t2_data$samples %in% common_samples, ]
            }
            results <- parallel::mclapply(name_combs, function(x) {
                  m1 <- t1_data[ ,x$V1]
                  m2 <- t2_data[ ,x$V2]
                  result <- rstatix::t_test(m1 ~ m2, paired=paired, alternative=type)
                  return(result)
              }, max.cores=4) %>% do.call(what=rbind.data.frame)
            return(results)
      })
        return(out)
  })


