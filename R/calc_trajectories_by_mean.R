
#' Function to get mean trajectories of metabolites and phenotypic traits 
#' @description function to extract mean trajectories
#' @param object An S4 object of class metime_analyser
#' @param which_data Dataset of interest
#' @param columns Other data that you want to see along with metabolites(column names from rowdata)
#' @return plot_data table that can be used to make the plotter object without metadata
#' @export
setGeneric("calc_trajectories_by_mean", function(object, which_data, columns) standardGeneric("calc_trajectories_by_mean"))
setMethod("calc_trajectories_by_mean", "metime_analyser", function(object, which_data, columns) {
      #Make sure the data is already scaled
      data <- object@list_of_data[[which_data]]
      data <- data[order(rownames(data)), ]
      if(!is.null(columns)) {
        rowdata <- object@list_of_row_data[[which_data]]
        rowdata <- rowdata[ ,columns]
        rowdata <- scale(rowdata, center=TRUE, scale=TRUE)
        rowdata <- rowdata[order(rownames(data)), ]
        final <- as.data.frame(cbind(data, rowdata))
      } else {
        final <- data
      }
      #final <- na.omit(final)
      object@list_of_data[[which_data]] <- final
      object <- mod_split_acc_to_time(object)
      data <- object@list_of_data[[which_data]]
      n_samples <- lapply(data, function(x) {
            return(dim(x)[1])
        })
      means <- lapply(data, function(x) {
            y <- apply(x, 2, mean)
            return(y)
        })
      means <- do.call(rbind, means)
      timepoints <- rownames(means)
      means <- as.data.frame(cbind(means, n_samples, timepoints))
      means$timepoints <- as.character(means$timepoints)
      means$n_samples <- as.character(means$n_samples)
      means <- means[order(as.numeric(rownames(means))), ]
      means <- as.data.frame(apply(means, 2, as.numeric))
      means <- reshape2::melt(means, id.vars=c("timepoints", "n_samples"))
      colnames(means) <- c("timepoints", "n_samples", "variables", "means")
      sds <- lapply(data, function(x) {
            y <- apply(x, 2, sd)
            return(y)
        })
      sds <- do.call(rbind, sds)
      #colnames(sds) <- paste(colnames(sds), "_sd", sep="")
      sds <- as.data.frame(cbind(sds, timepoints))
      sds <- sds[order(as.numeric(rownames(sds))), ]
      sds <- as.data.frame(apply(sds, 2, as.numeric))
      sds <-  reshape2::melt(sds, id.vars=c("timepoints"))
      colnames(sds) <- c("timepoints", "variables", "sd")
      sd <- sds$sd
      vars <- lapply(data, function(x) {
            y <- apply(x, 2, var)
            return(y)
        })
      vars <- do.call(rbind, vars)
      #colnames(vars) <- paste(colnames(vars), "_var", sep="")
      vars <- as.data.frame(cbind(vars, timepoints))
      vars <- vars[order(as.numeric(rownames(vars))), ] 
      vars <- as.data.frame(apply(vars, 2, as.numeric))
      vars <- reshape2::melt(vars, id.vars=c("timepoints"))
      colnames(vars) <- c("timepoints", "variables", "var")
      variance <- vars$var
      out <- as.data.frame(cbind(means, sd, variance))
      return(out)  
  })
 
