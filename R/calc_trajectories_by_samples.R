#add standard deviation and variance
#zscore the column of interest and the raw data to line plot the trajectory 

#calc_trajectories_by_samples
setGeneric("calc_trajectories_by_samples", function(object, which_data, columns) standardGeneric("calc_trajectories_by_samples"))
setMethod("calc_trajectories_by_samples", "metime_analyser", function(object, which_data, columns) {
      #Make sure the data is already scaled
      data <- object@list_of_data[[which_data]]
      rowdata <- object@list_of_row_data[[which_data]]
      rowdata <- rowdata[ ,columns]
      rowdata <- scale(rowdata, center=TRUE, scale=TRUE)
      data <- data[order(rownames(data)), ]
      rowdata <- rowdata[order(rownames(data)), ]
      final <- as.data.frame(cbind(data, rowdata))
      #final <- na.omit(final)
      data <- data[order(rownames(data)), ]
      samples <- as.character(unlist(lapply(strsplit(rownames(final), split="_"), function(x) return(x[1]))))
      indList <- split(seq_along(samples), samples)
      #sample wise separation of the data
      list_of_samples_data <- lapply(indList, function(x) {
          return(data[x,])
      })
      list_of_samples_data <- lapply(list_of_samples_data, function(x) {
          timepoints <- as.numeric(unlist(lapply(strsplit(rownames(x), split="_"), function(x) return(x[2]))))
          samples <- as.character(unlist(lapply(strsplit(rownames(x), split="_"), function(x) return(x[1]))))
          x <- as.data.frame(cbind(x, timepoints, samples))
          x <- x[order(as.numeric(x$timepoints)), ]
          return(x)
      })
      plot_data <- do.call(rbind, list_of_samples_data)
      plot_data <- reshape2::melt(plot_data, id.vars=c("samples", "timepoints"))
      return(plot_data)
  }) 

