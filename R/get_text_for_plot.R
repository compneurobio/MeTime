

#' Function to Obtain textual information for visualization in interactive plots
#' @description a standard function to be applied on data matrices or dataframes with the colnames of interest such that the information from
#' columns is visualized in the interactive plot
#' @examples
#' # text = get_text(data=data.frame, colnames=c("names","of","columns", "of", "interest"))
#' @param data a dataframe with plotting data along with other variables for visualization
#' @param colnames a character vector with the names of the variables that you want to see on the plot
#' @return a vector with strings that can be parsed into plot_ly text.
#' @export
get_text_for_plot <- function(data, colnames) {	
		out <- c()
		count <- 1
		text <- c()
		for(l in 1:length(rownames(data))) {
			for(m in 1:length(colnames)) {
				if(m==1) {
						text <- paste("<br /> ", colnames[m], " : ", data[l, colnames[m]], sep="")
					} else {
						text <- paste(text, "<br /> ", colnames[m], " : ", data[l, colnames[m]], sep="")
					}
				}
				out[count] <- text
				count <- count + 1
			} 

		return(out)
	} 

