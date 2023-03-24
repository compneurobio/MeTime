
#' Function to Convert S4 object of class metime_analyser to an S3 object(list) with same architecture
#' @description Converter function to be applied onto metime_analyse object to convert into a standard list of S3 type.
#' For example list_of_data can now be accessed with '$' instead of '@'
#' @examples # convert S4 object to a list
#' s3_list <- mod_convert_s4_to_s3(object=metime_analyser_object)
#' @param object An object of class metime_analyser
#' @return A list with the same data and structure as a metime_analyser 
#' @export
setGeneric("mod_convert_s4_to_s3", function(object) standardGeneric("mod_convert_s4_to_s3"))

setMethod("mod_convert_s4_to_s3", "metime_analyser", function(object) {
		#will add based on analysis - make it module wise or open for suggestions
		out <- list(list_of_data=object@list_of_data, list_of_col_data=object@list_of_col_data, 
			list_of_row_data=object@list_of_row_data, annotations=object@annotations,
			results=object@results)
		return(out)
	})


