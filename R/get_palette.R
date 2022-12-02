#' Get a palette of "n" distinct colorblind friendly colors 
#' @description Function to get a palette of distinct colorblind friendly colors, the distinctiveness is determined by the difference in their hue values.
#' @examples
#' # colors=get_palette(n=10)
#' @param n number of colors wanted in the palette
#' @return a color palette vector with colors in the form of hex codes
#' @export 

get_palette <- function(n) {
	#loading the package to get colors 
	#require(RColorBrewer)
	#extracting all the colorblind friendly colors
	colors <- RColorBrewer::brewer.pal.info
	colors <- colors[colors$colorblind == TRUE, ]
	col_vec <- unlist(mapply(RColorBrewer::brewer.pal, colors$maxcolors, rownames(colors)))
	col_vec <- unique(col_vec)
	#get distinct colors by converting the hex to rgb and then to HSL values
	#We then order the colors based on hue and the differences between each are ranked
	#we then retain the 'n' colors needed
	hue_df <- data.frame(color=col_vec, hue=rgb2hsv(col2rgb(col_vec))[1,])#here first row is hue values
	hue_df <- hue_df[order(hue_df$hue),]
  	hue_df <- rbind(hue_df, hue_df[1,])
  	hue_df$hue[nrow(hue_df)] <- hue_df$hue[nrow(hue_df)] + 1
  	hue_df$diff <- c(1, diff(hue_df$hue))
  	hue_df <- hue_df[order(-hue_df$diff),]
  	hue_df <- hue_df[-1, ]
  	out <- hue_df$color[1:n]
  	return(out)
}

