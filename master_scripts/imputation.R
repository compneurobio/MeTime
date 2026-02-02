require(MeTime)

#loading the imputed analyser object
load("adni_q500_data")
# this data is already imputed so we will randomly introduce missingness
# to showcase this example
df <- get_data(adni_q500_data, which_data="biocrates_data")
sample_names <- rownames(df)
df <- as.data.frame(lapply(df, function(cc) cc[ sample(c(TRUE, NA), prob = c(0.85, 0.15), size = length(cc), replace = TRUE) ]))
rownames(df) <- sample_names
adni_q500_data@list_of_data$nmr_data <- df

# checking the missingness and plotting for random 100 metabolites
adni_q500_data %>% plot_missingness(which_data="biocrates_data", max_features=100)

# checking thresholds
adni_q500_data <- adni_q500_data %>% 
                    mod_filter_features_by_missingness(which_data="biocrates_data", threshold=0.3) %>%
                    mod_filter_samples_by_missingness(which_data="biocrates_data", threshold=0.3) %>%
                    mod_filter(which_data="biocrates_data", type="row_data", BIFAST==1) %>%
                    #mod_trans_log(which_data="nmr_data", base=2) %>% currently commented out as data is already log transformed
                    mod_impute(which_data = which_data,
                     thresh = 0.3,
                     path_to_diagnostics = getwd(), # change this to directory of interest
                     method="rf")
