# Package index

## All functions

- [`add_data()`](https://compneurobio.github.io/MeTime/reference/add_data.md)
  : Add a data frame to analyzer object.
- [`add_dataset()`](https://compneurobio.github.io/MeTime/reference/add_dataset.md)
  : This function appends an object of class metime_analyser with a new
  dataset.
- [`add_distribution_vars_to_rows()`](https://compneurobio.github.io/MeTime/reference/add_distribution_vars_to_rows.md)
  : Add phenotypic measurements that are not added to the row_data of
  the dataset
- [`add_function_info()`](https://compneurobio.github.io/MeTime/reference/add_function_info.md)
  : Add information of function added to the data
- [`add_node_features()`](https://compneurobio.github.io/MeTime/reference/add_node_features.md)
  : Function to add features to visnetwork plot from another plotter
  object
- [`add_plot()`](https://compneurobio.github.io/MeTime/reference/add_plot.md)
  : Generates a report from an metime_analyser object
- [`add_result()`](https://compneurobio.github.io/MeTime/reference/add_result.md)
  : Add a result element to a metime_analyser
- [`add_screening_vars()`](https://compneurobio.github.io/MeTime/reference/add_screening_vars.md)
  : Add values measured at a screening time for samples to be added to
  all time points
- [`calc_baseline_change()`](https://compneurobio.github.io/MeTime/reference/calc_baseline_change.md)
  : Calculate baseline change per subject
- [`calc_clusters_wgcna()`](https://compneurobio.github.io/MeTime/reference/calc_clusters_wgcna.md)
  : Calculate cluster assignment from WGCNA
- [`calc_col_stats()`](https://compneurobio.github.io/MeTime/reference/calc_col_stats.md)
  : Check normality of metabolites and append results to col_data
- [`calc_colinearity_features()`](https://compneurobio.github.io/MeTime/reference/calc_colinearity_features.md)
  : Function to calculate colinearity of features
- [`calc_conservation_metabolite()`](https://compneurobio.github.io/MeTime/reference/calc_conservation_metabolite.md)
  : Calculation function to analyze the metabolite conservation index
- [`calc_conservation_metabotype()`](https://compneurobio.github.io/MeTime/reference/calc_conservation_metabotype.md)
  : Function to calculate metabotype conservation index
- [`calc_correlation_features()`](https://compneurobio.github.io/MeTime/reference/calc_correlation_features.md)
  : Function to calculate correlation
- [`calc_dimensionality_reduction_metabs()`](https://compneurobio.github.io/MeTime/reference/calc_dimensionality_reduction_metabs.md)
  : Function to calculate dimensionality reduction methods such as tsne,
  umap and pca.
- [`calc_dimensionality_reduction_samples()`](https://compneurobio.github.io/MeTime/reference/calc_dimensionality_reduction_samples.md)
  : Function to calculate dimensionality reduction methods such as tsne,
  umap and pca.
- [`calc_distance_samples()`](https://compneurobio.github.io/MeTime/reference/calc_distance_samples.md)
  : Function to calculate dissimilarity using distance measures
- [`calc_distribution_metabs()`](https://compneurobio.github.io/MeTime/reference/calc_distribution_metabs.md)
  : Function for Plotting distributions of phenotypic variables
- [`calc_distribution_samples()`](https://compneurobio.github.io/MeTime/reference/calc_distribution_samples.md)
  : Function for Plotting distributions of phenotypic variables
- [`calc_featureselection_boruta()`](https://compneurobio.github.io/MeTime/reference/calc_featureselection_boruta.md)
  : Function to calculate dependent variables
- [`calc_gamm()`](https://compneurobio.github.io/MeTime/reference/calc_gamm.md)
  : Calculation of generalized additive mixed models (GAMMs)
- [`calc_ggm_genenet()`](https://compneurobio.github.io/MeTime/reference/calc_ggm_genenet.md)
  : An automated fucntion to calculate GGM from genenet crosssectional
  version
- [`calc_ggm_multibipartite_lasso()`](https://compneurobio.github.io/MeTime/reference/calc_ggm_multibipartite_lasso.md)
  : An automated fucntion to calculate GGM from multibipartite lasso
  approach
- [`calc_lm()`](https://compneurobio.github.io/MeTime/reference/calc_lm.md)
  : Cossectional linear models (per time point)
- [`calc_lmm()`](https://compneurobio.github.io/MeTime/reference/calc_lmm.md)
  : Calculation of linear mixed models
- [`calc_temporal_network()`](https://compneurobio.github.io/MeTime/reference/calc_temporal_network.md)
  : An automated function to caluclate temporal network with lagged
  model
- [`calc_time_trend()`](https://compneurobio.github.io/MeTime/reference/calc_time_trend.md)
  : Calculate feature time trends
- [`check_col_normality()`](https://compneurobio.github.io/MeTime/reference/check_col_normality.md)
  : Function to check for col_normality data whether it is added or not.
- [`check_scaling_and_transformation()`](https://compneurobio.github.io/MeTime/reference/check_scaling_and_transformation.md)
  : Function to check if the data is already scaled or log transformed
- [`get_class_info_from_edges()`](https://compneurobio.github.io/MeTime/reference/get_class_info_from_edges.md)
  : Get summary on class edges
- [`get_coldata()`](https://compneurobio.github.io/MeTime/reference/get_coldata.md)
  : Get col_data from a S4 object of the class "metime_analyser"
- [`get_common_samples_at_timepoints()`](https://compneurobio.github.io/MeTime/reference/get_common_samples_at_timepoints.md)
  : Get common samples at multiple timepoints chosen
- [`get_convert_to_se()`](https://compneurobio.github.io/MeTime/reference/get_convert_to_se.md)
  : Convert analyser datasets to SummarizedExperiment objects
- [`get_data()`](https://compneurobio.github.io/MeTime/reference/get_data.md)
  : Get data from a S4 object of the class "metime_analyser"
- [`get_files_and_names()`](https://compneurobio.github.io/MeTime/reference/get_files_and_names.md)
  : Pack a dataset (data, row_data, col_data) into a single object of
  class "metime_analyser".
- [`get_li_thresh()`](https://compneurobio.github.io/MeTime/reference/get_li_thresh.md)
  : Calculate significance threshold after correcting for independent
  tests.
- [`get_make_analyser_object()`](https://compneurobio.github.io/MeTime/reference/get_make_analyser_object.md)
  : Pack all the data into a single object of class "metime_analyser"
- [`get_make_results()`](https://compneurobio.github.io/MeTime/reference/get_make_results.md)
  : Get results list for S4 object of class "metime_analyser" object
- [`get_metadata_for_columns()`](https://compneurobio.github.io/MeTime/reference/get_metadata_for_columns.md)
  : Get metadata for columns(in most cases for metabolites)
- [`get_metadata_for_rows()`](https://compneurobio.github.io/MeTime/reference/get_metadata_for_rows.md)
  : Get metadata for rows(in most cases for samples)
- [`get_parameters_of_results()`](https://compneurobio.github.io/MeTime/reference/get_parameters_of_results.md)
  : Get parameters of the results
- [`get_results()`](https://compneurobio.github.io/MeTime/reference/get_results.md)
  : Get results element by index
- [`get_rowdata()`](https://compneurobio.github.io/MeTime/reference/get_rowdata.md)
  : Get row_data from a S4 object of the class "metime_analyser"
- [`get_samples_and_timepoints()`](https://compneurobio.github.io/MeTime/reference/get_samples_and_timepoints.md)
  : Summarize number of time points and the total number of samples
  available at that point
- [`get_stratified_data()`](https://compneurobio.github.io/MeTime/reference/get_stratified_data.md)
  : Function to stratify before calculation
- [`meta_conservation()`](https://compneurobio.github.io/MeTime/reference/meta_conservation.md)
  : Meta comparison for conservation index results
- [`meta_feature_overlap()`](https://compneurobio.github.io/MeTime/reference/meta_feature_overlap.md)
  : Meta comparison for feature overlap
- [`meta_matrix_similarity()`](https://compneurobio.github.io/MeTime/reference/meta_matrix_similarity.md)
  : Meta comparison for matrix similarity
- [`meta_network_overlap()`](https://compneurobio.github.io/MeTime/reference/meta_network_overlap.md)
  : Meta comparison for network overlap
- [`meta_regression()`](https://compneurobio.github.io/MeTime/reference/meta_regression.md)
  : Meta comparison for regression outputs
- [`metime_analyser-class`](https://compneurobio.github.io/MeTime/reference/metime_analyser.md)
  : Constructor to generate an object of class metime_analyser. contains
  slots - list_of_data: For the list of all data matrices. -
  list_of_col_data: list of all the col data files in the same order. -
  list_of_row_data: list of all the row data files in the same order. -
  annotations: list with phenotype and medication. Each of which is
  character that represents the name of the aforementioned dataset
  types.
- [`mod_convert_s4_to_s3()`](https://compneurobio.github.io/MeTime/reference/mod_convert_s4_to_s3.md)
  : Function to Convert S4 object of class metime_analyser to an S3
  object(list) with same architecture
- [`mod_extract_common_samples()`](https://compneurobio.github.io/MeTime/reference/mod_extract_common_samples.md)
  : Function to extract common samples across all datasets and store
  them only
- [`mod_filter()`](https://compneurobio.github.io/MeTime/reference/mod_filter.md)
  : Modification function to filter columns in data, row_data or
  col_data
- [`mod_filter_features_by_missingness()`](https://compneurobio.github.io/MeTime/reference/mod_filter_features_by_missingness.md)
  : Filter features by missingness
- [`mod_filter_features_by_variance()`](https://compneurobio.github.io/MeTime/reference/mod_filter_features_by_variance.md)
  : Filter features by variance
- [`mod_filter_samples_by_missingness()`](https://compneurobio.github.io/MeTime/reference/mod_filter_samples_by_missingness.md)
  : Filter samples by missingness
- [`mod_filter_subjects()`](https://compneurobio.github.io/MeTime/reference/mod_filter_subjects.md)
  : Filter subjects from datasets
- [`mod_filter_timepoints()`](https://compneurobio.github.io/MeTime/reference/mod_filter_timepoints.md)
  : Filter timepoints from datasets
- [`mod_generate_plots()`](https://compneurobio.github.io/MeTime/reference/mod_generate_plots.md)
  : Function to update plots post calculations
- [`mod_impute()`](https://compneurobio.github.io/MeTime/reference/mod_impute.md)
  : Imputation of missing values
- [`mod_merge_data()`](https://compneurobio.github.io/MeTime/reference/mod_merge_data.md)
  : Modification (mod) function to merge two sets of data or partial
  data merging for any analysis
- [`mod_merge_metime_analysers()`](https://compneurobio.github.io/MeTime/reference/mod_merge_metime_analysers.md)
  : Function to merge one or more metime_analyser objects
- [`mod_merge_results()`](https://compneurobio.github.io/MeTime/reference/mod_merge_results.md)
  : Function to merge two sets of results
- [`mod_merge_row_data_and_data()`](https://compneurobio.github.io/MeTime/reference/mod_merge_row_data_and_data.md)
  : Modification (mod) function to merge data and row_data (sample info)
  partially or completely.
- [`mod_mutate()`](https://compneurobio.github.io/MeTime/reference/mod_mutate.md)
  : Function to mutate columns in row_data or col_data
- [`mod_remove_duplicates()`](https://compneurobio.github.io/MeTime/reference/mod_remove_duplicates.md)
  : Function to remove duplicates
- [`mod_remove_nas()`](https://compneurobio.github.io/MeTime/reference/mod_remove_nas.md)
  : Function to remove NA's from data matrices
- [`mod_rename()`](https://compneurobio.github.io/MeTime/reference/mod_rename.md)
  : Function to rename columns in the datasets and results of analyser
  object
- [`mod_select()`](https://compneurobio.github.io/MeTime/reference/mod_select.md)
  : Modification function to select columns in data, row_data or
  col_data
- [`mod_trans_eigendata()`](https://compneurobio.github.io/MeTime/reference/mod_trans_eigendata.md)
  : Function to add the clusters obtained from wgcna
- [`mod_trans_foldchange()`](https://compneurobio.github.io/MeTime/reference/mod_trans_foldchange.md)
  : Function to calculate fold change
- [`mod_trans_log()`](https://compneurobio.github.io/MeTime/reference/mod_trans_log.md)
  : Function to apply log transformation
- [`mod_trans_zscore()`](https://compneurobio.github.io/MeTime/reference/mod_trans_zscore.md)
  : Function to scale the data
- [`plot(`*`<metime_analyser>`*`,`*`<ANY>`*`)`](https://compneurobio.github.io/MeTime/reference/plot-metime_analyser-ANY-method.md)
  : Setting a plotting method for the metime_analyser class
- [`plot_missingness()`](https://compneurobio.github.io/MeTime/reference/plot_missingness.md)
  : Plot missingness heatmap for a dataset
- [`rm_data()`](https://compneurobio.github.io/MeTime/reference/rm_data.md)
  : Function to remove datasets or results from the metime analyser
  object
- [`show(`*`<metime_analyser>`*`)`](https://compneurobio.github.io/MeTime/reference/show-metime_analyser-method.md)
  : Setting new print definition for the metime_analyser object
- [`summarize_dataset()`](https://compneurobio.github.io/MeTime/reference/summarize_dataset.md)
  : Summarize datasets in a metime_analyser object
- [`validate_metime_analyser()`](https://compneurobio.github.io/MeTime/reference/validate_metime_analyser.md)
  : Validate a metime_analyser object
- [`validity()`](https://compneurobio.github.io/MeTime/reference/validity.md)
  : Validity function to check if the object is valid or not
- [`write_data()`](https://compneurobio.github.io/MeTime/reference/write_data.md)
  : Function to extract analyser object data into a csv
- [`write_report()`](https://compneurobio.github.io/MeTime/reference/write_report.md)
  : Generates a report from an metime_analyser object
- [`write_results()`](https://compneurobio.github.io/MeTime/reference/write_results.md)
  : Function to save plot_data from an analyser_objec into an xlsx file
