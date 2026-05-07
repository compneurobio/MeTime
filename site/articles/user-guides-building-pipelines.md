# Building Pipelines

## Building Pipelines

> One-line summary: chain modular functions into reproducible analysis
> pipelines.

### Pipeline mental model

1.  Build object
2.  Modify data (`mod_*`)
3.  Calculate (`calc_*`)
4.  Plot/export
5.  Meta-analysis (`meta_*`, optional)

For most of the calculation external data is not needed but for
regression analysis an additional dataset has to be appended to
col_data. See below for example in humet and in case studies

### Starter pipeline template ()

``` r
# load package
library(MeTime)
# build analyser object
obj <- get_make_analyser_object(...)
# perform calculation and generate plots
obj <- mod_<step1>(...) %>%
        mod_<step2>(...) %>%
        calc_<analysis>(which_data = "<dataset>", name = "<result_1>") %>% 
        mod_generate_plots(...)
# write report for the results
obj %>% write_report(...)
```

### Pipelines for each type of question

Here are a few examples for building pipelines see case studies example
for more details

1.  Distribution, distances and correlation

Here is an example workflow covering distributions, correlations, and
distances with HuMet data

``` r

data("humet_object")
# data distributions
which_data <- "humet_subset_data"
distributions <- humet_object %>%
        mod_trans_zscore(which_data=which_data) %>%
        calc_distribution_samples(which_data=which_data, 
            cols=c("Factor.Value.Challenge.", "Factor.Value.Identifier."), stratifications=NULL, name="samples_distribution") %>%
        mod_generate_plots(type="distribution_samples", .interactive=TRUE) %>%
        calc_distribution_metabs(which_data=which_data, 
            cols=c("super_pathway", "sub_pathway"),name="metabs_distribution") %>%
        mod_generate_plots(type="distribution_metabs", .interactive=TRUE)

distances <- humet_object %>% 
    mod_trans_zscore(which_data=which_data) %>%
    calc_col_stats(which_data=which_data, cols_for_meta=NULL) %>%
    calc_correlation_features(which_data=which_data, method="spearman", 
                    stratifications=NULL) %>%
    mod_generate_plots(type="pairwise_correlation", .interactive=T) %>% 
    calc_distance_samples(which_data=which_data, method="euclidean", 
                    stratifications=NULL) %>%
    mod_generate_plots(type="pairwise_distance", .interactive=TRUE)
```

2.  Dimensionality reduction

Example already shown in quickstart. The approach is the same for tSNE
and UMAP

3.  Conservation index

Example already shown in data_preparation

4.  Imputation and feature selection

``` r
which_data <- "humet_data"
imputed_object <- humet_object %>%
        # removing samples with more than >30% missingness
        mod_filter_samples_by_missingness(threshold=0.3, which_data=which_data) %>%
        # removing features with more than >30% missingness
        mod_filter_features_by_missingness(which_data=which_data, threshold=0.3) %>%
        # imputing data
        mod_impute(method="rf", which_data=which_data)

# for feature selection example, creating a dummy medication dataset for showcase
# and also using the subset data for ease
which_data <- "humet_subset_data"
data <- get_data(humet_object, which_data=which_data)
n_samples <- nrow(data)
n_meds <- 50
# Use sample IDs from your existing dataset
sample_ids <- rownames(data)
# Medication IDs
med_ids <- sprintf("MED_%03d", 1:n_meds)
# Dummy medication intake matrix
med_data <- matrix(
  rbinom(n_samples * n_meds, size = 1, prob = 0.2),
  nrow = n_samples,
  ncol = n_meds,
  dimnames = list(sample_ids, med_ids)
) %>% as.data.frame()
disease_categories <- c(
  "Hypertension",
  "Diabetes",
  "Hyperlipidemia",
  "Depression",
  "Pain",
  "Cardiovascular disease",
  "Inflammation",
  "Asthma",
  "Alzheimer's disease",
  "Parkinson's disease",
  "Anxiety",
  "Infection"
)

med_metadata <- data.frame(
  ids = med_ids,
  pathway_or_disease = sample(disease_categories, n_meds, replace = TRUE),
  stringsAsFactors = FALSE
)

feature_selection <- humet_object %>%
    add_dataset(data=med_data, col_data=med_metadata, row_data=get_rowdata(humet_object, which_data="humet_subset_data"), name="medication_data") %>%
    mod_trans_zscore(which_data="humet_subset_data") %>%
# this part is computationally expensive so beware
#  calc_featureselection_boruta(which_x="medication_data", # here we use medication but you can use any other metabolomics dataset
#                                       which_y="humet_subset_data",
#                                       verbose=T,
#                                       name=paste0(which_data, "_selected_full"),
#                                       cols_for_meta_x=NULL,
#                                       cols_for_meta_y=NULL,
#                                       save_per_run=T,
#                                       num_cores=12)
```

5.  Calculating clusters and eigendata

``` r

which_data <- "humet_subset_data"
eigendata <- humet_object %>%
  mod_trans_zscore(which_data=which_data) %>% 
  mod_trans_eigendata(which_data = which_data, append=T, 
                      cols_for_meta=list(humet_subset_data=c(id="id", sub_pathway="sub_pathway")),
                      baseline="t1",
                      minClusterSize=2,
                      deepSplit=3,
                      pamRespectsDendro=TRUE) 

# you can also store it as a result if you don't want to use the eigendata for further analysis

eigen_result <- humet_object %>%
    mod_trans_zscore(which_data=which_data) %>% 
    calc_clusters_wgcna(which_data = which_data, baseline="t1",
                      cols_for_meta=list(humet_subset_data=c(id="id", sub_pathway="sub_pathway")),
                      minClusterSize=2,
                      deepSplit=3,
                      pamRespectsDendro=TRUE)
```

6.  Network analysis

Example for partial correlation based network

``` r

which_data <- "humet_subset_data"
networks <- humet_object %>%
  mod_merge_row_data_and_data(which_data=which_data, 
        cols_list=list(data=NULL, row_data=c("Factor.Value.Challenge.")),
        name="ggm_data", append=T) %>%
  calc_ggm_genenet(which_data = "ggm_data", threshold = "li", all=FALSE, 
                     cols_for_meta = list(humet_subset_data=c(sub_pathway="sub_pathway")),
                     covariates = c("Factor.Value.Challenge."),
                     stratifications = list(time=c("t1", "t2", "t3", "t4", "t5")),
                     name="genenet_ggm_results") %>%
    mod_generate_plots(type="network")
```

7.  Regression analysis

Several different kinds of mixed models can be implemented using this
package we show you the simplest linear mixed model example here. For
more complicated analysis please look at the case studies.

``` r

data("humet_object")
which_data <- "humet_subset_data"
col_data <- get_coldata(humet_object, which_data=which_data)

lmm_data <- col_data %>% select(id) %>% mutate(cov=NA, type="met", interaction=NA, random=NA) 
traits <- data.frame(id="Factor.Value.Challenge.", cov=NA, type="trait", random="subject", interaction=NA)
lmm_data <- plyr::rbind.fill(lmm_data, traits)

# the goal is to have a dataframe with id, covs, type, interaction and random in case of linear mixed models and generalized additive models whereas for linear model interaction and random are not needed.
# For cases where you need more than one covariate or random effect or interaction effect just add the variable names with ### separating them for example: ###Age###BMI. Same logic can be applied to cov column as well
# In this example we created a dataframe but we recommend to use an excel sheet
# An example excel sheet can be found in the case-studies folder

linear_mixed_model <- humet_object %>%
    mod_merge_row_data_and_data(which_data=which_data, 
        cols_list=list(data=NULL, row_data=c("Factor.Value.Challenge.")),
        name="lmm_data", append=T) %>%
    add_data(which_data="lmm_data", type="col_data",x=lmm_data, id="id") %>%
    calc_lmm(which_data="lmm_data", 
             name="simple_mixed_model", 
             stratifications=NULL, 
            cols_for_meta=list(c(sub_pathway="sub_pathway")), # change column name as per need
             num_cores=12)
```

If the user is using a windows operating system, PSOCK might create
overthreading and to prevent that set these environment variables before
hand in R for a faster performance. In case of linux based operating
systems this is not needed.

``` r

Sys.setenv(OMP_NUM_THREADS=1)
Sys.setenv(MKL_NUM_THREADS=1)
Sys.setenv(OPENBLAS_NUM_THREADS=1)
```
