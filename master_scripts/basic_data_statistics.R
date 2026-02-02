require(MeTime)

#loading the imputed analyser object
load("adni_q500_data.rda")
which_data <- "biocrates_data"
distrubutions <- adni_q500_data %>%
	add_distribution_vars_to_rows(screening_vars=NULL, 
                                distribution_vars=c("APOEGrp", "PTGENDER", "Age", "BMI", "PTEDUCAT", "DXGrp_long"), 
                                which_data=which_data) %>%
 	calc_distribution_samples(which_data=which_data, cols=c("APOEGrp",  "DXGrp_long", "Age", "BMI", "PTEDUCAT"),
                             stratifications=NULL, name="samples_distribution") %>%
   	mod_generate_plots(type="distribution_samples", .interactive=TRUE) %>%
 	calc_distribution_metabs(which_data=which_data, cols=c("Class"), name="metabs_distribution") %>%
   	mod_generate_plots(type="distribution_metabs", .interactive=TRUE)

distrubutions %>% 
   	write_report(file="distributions.html", title="ADNI Q500 data variables information")


distances <- adni_q500_data %>% 
			calc_col_stats(which_data=which_data, cols_for_meta=NULL) %>%
			calc_correlation_features(which_data=which_data, method="spearman", 
					stratifications=NULL) %>%
			mod_generate_plots(type="pairwise_correlation", .interactive=T) %>% 
			calc_distance_samples(which_data=which_data, method="euclidean", 
					stratifications=NULL) %>%
			mod_generate_plots(type="pairwise_distance", .interactive=TRUE)


# Change the plotting to non-interactive style to generate the markdowns 
# otherwise the size of the plots is too big to generate plots
#distances %>%
#	write_report(file="distances.html", title="ADNI Q500 data correlations")

	