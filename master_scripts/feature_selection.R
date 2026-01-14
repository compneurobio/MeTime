require(MeTime)

# Example feature selection workflow with optional medication covariates.

# Loading the imputed analyser object.
load("adni_nmr_data")
# Dataset of interest.
dataset_name <- "nmr_data"
which_data <- dataset_name
medication_name <- "medication_data"
exclude_meds <-  c("N06DX","N06DA") # exclude anti-dementia drugs

# Prepare data for feature selection.
feature_selection <- adni_nmr_data %>% 
				 	add_distribution_vars_to_rows(screening_vars=NULL, 
							distribution_vars=c("ADNI_MEM", "ADNI_LAN", "ADNI_EF", 
								"APOEGrp",  "DXGrp_longi", "PTGENDER", "Age", "BMI"), 
							which_data=which_data) %>% 
				 	mod_mutate(which_data=which_data, type="row_data", 
				 		PTGENDER=as.factor(PTGENDER), APOEGrp=as.factor(APOEGrp)) %>% 
          mod_rename(which_data = which_data, type="col_data", super_pathway="Group", sub_pathway="Subgroup") %>% 
          mod_select(which_data= which_data, type="data", !contains("pct|PCT")) 

# In this example we showcase selection of medication for later using them as covariates.
# DO NOT RUN this locally; it is commented out.
# Run feature selection on the three main timepoints:
# selected_full <- feature_selection %>% 
#   MeTime::mod_filter(time %in% c("t0","t12","t24"), which_data=dataset_name, type="row_data") %>% # filter data to main time points
#   MeTime::mod_filter(time %in% c("0","12","24"), which_data=medication_name, type="row_data") %>% # filter med_data to main time points
#   MeTime::mod_filter(
#     id %in% setdiff(feature_selection@list_of_col_data[[medication_name]]$id, exclude_meds),
#     which_data=medication_name,
#     type="col_data"
#   ) %>% # filter any medications that should be excluded
#   MeTime::calc_featureselection_boruta(
#     which_x=medication_name, # here we use medication but you can use any other dataset
#     which_y=dataset_name,
#     verbose=TRUE,
#     name=paste0(dataset_name, "_selected_full"),
#     cols_for_meta_x=NULL,
#     cols_for_meta_y=NULL,
#     save_per_run=TRUE,
#     num_cores=31 # change the cores as per need
#   )
# saveRDS(selected_full, file=paste0(path_out, "metabolon_medsel_all.rds"))



