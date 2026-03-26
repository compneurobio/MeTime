# Function Families

> One-line summary: use function prefixes to navigate the API quickly.

## `add_*`
- Purpose: append new elements.
- Typical usage: add datasets, plots, metadata.
- Examples: `add_data`, `add_result`, `add_plot`.

## `mod_*`
- Purpose: mutate/filter/transform data.
- Typical usage: preprocessing and harmonization.
- Examples: `mod_filter`, `mod_impute`, `mod_trans_log`.
- Wrappers: For easier data manipulation we have wrappers for mutate, select, rename and fiter from the dplyr package
- QC helpers: There are additional mod_ functions that can help user to remove features/samples based on missingness or on variance (Example: `mod_filter_samples_by_missingness`)

## `calc_*`
- Purpose: perform analysis methods.
- Typical usage: one result object per named analysis.
- These set of analysis are currently available, a more detailed description can be found in choosing methods:
    - Distributions, distances and correlations
    - Feature selection
    - Imputation
    - Dimensionality reduction
    - Eigendata calculation
    - Conservation index analysis
    - Regression (Linear models, Linear mixed models and Generalized additive mixed models)
    - Data-driven networks 
    

## `meta_*`
- Purpose: compare or aggregate comparable results.
- Functions include:
    - meta_regression (compare results of similar models)
    - meta_conservation (compare conservation index analysis across two results)
    - meta_feature_overlap (compare feature selection analysis across two results)
    - meta_network_overlap (compare two networks)
    - meta_matrix_similarity (compare two results from calc_distance_samples or calc_correlation_metabolites)

## `get_*`
- Purpose: inspect and extract data/results.

## `write_*`
- Purpose: export results and reports.
