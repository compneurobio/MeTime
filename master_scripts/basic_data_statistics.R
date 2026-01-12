require(MeTime)

# Example workflow for basic data statistics in MeTime.
# This script assumes that "adni_nmr_data" has already been imputed and saved.

# Loading the imputed analyser object.
load("adni_nmr_data")
which_data <- "nmr_data"

# 1) Distribution summaries for samples and metabolites.
adni_nmr_data <- adni_nmr_data %>%
	add_distribution_vars_to_rows(screening_vars=NULL, 
                                distribution_vars=c("APOEGrp", "PTGENDER", "Age", "BMI", "PTEDUCAT", "DXGrp_long"), 
                                which_data=which_data) %>%
 	calc_distribution_samples(which_data=which_data, cols=c("ADNI_MEM", "ADNI_LAN", "ADNI_EF", "APOEGrp",  "DXGrp_longi", "Age", "BMI", "PTEDUCAT"),
                             stratifications=NULL, name="samples_distribution") %>%
   	mod_generate_plots(type="distribution") %>%
 	calc_distribution_metabs(which_data=which_data, cols=c("Group", "Subgroup"), name="metabs_distribution") %>%
   	mod_generate_plots(type="distribution")

adni_nmr_data %>% 
   	write_report(file="distributions.html", title="ADNI NMR data variables information")


# 2) Column statistics and correlation summaries.
adni_nmr_data <- adni_nmr_data %>% 
			calc_col_stats(which_data=which_data) %>%
			calc_correlation_pairwise(which_data=which_data, method="spearman", stratifications=NULL, 
			                          cols_for_meta=list(nmr_data=c(id="id", sub_pathway="Group"))) %>%
			mod_generate_plots(type="correlation")

adni_nmr_data %>% 
   	write_report(file="correlations.html", title="ADNI NMR correlation summaries")
