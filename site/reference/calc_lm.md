# Cossectional linear models (per time point)

Fit multiple linear models (lm) onto one dataset. Covariates can be
added using the covariate column of the col_data (multiple covariates
can be added, separated by '###') See examples for more information.

## Usage

``` r
calc_lm(
  object,
  which_data,
  stratifications = NULL,
  cols_for_meta = NULL,
  threshold = c("none", "nominal", "li", "fdr", "bonferroni"),
  verbose = T,
  name = "regression_lm_1",
  timepoint = NULL,
  num_cores = NULL,
  cluster_profile = c("auto", "local", "hpc"),
  hpc_libpaths = NULL
)
```

## Arguments

- object:

  an S4 object of class metime_analyser

- which_data:

  a character defining the name of the dataset to be used.

- stratifications:

  list to stratify data into a subset. Usage list(name=value). Default
  set to NULL, thereby not performing any type of stratification.

- cols_for_meta:

  a character vector to define column names that are to be used for
  plotting purposes. Default set to NULL, therby not adding columns as
  metadata. If you want automated facet wrapping option then set your
  new_columns as "facet_your_name"

- threshold:

  a character vector to define the type of threshold for significant
  interactions. Default set to all availabe thresholds:
  c("none","nominal","li","fdr","bonferroni"). allowed inputs are "li",
  "FDR", "bonferroni" and "nominal"(cutoff p=0.05, set as Default)

- verbose:

  a logical on whether to print the calculation progress. Default set to
  FALSE.

- name:

  a character vector to define the index within the results. Should be
  equal to length of which_data. Default set to regression_lm_1.

- timepoint:

  time input for cross-sectional model should be the same as the value
  in time column in row_data.

- num_cores:

  numeric input to define the number of cores that you want to use for
  parallel computing. Default is set to NULL which is
  parallel::detectCores() -1.

- cluster_profile:

  Character string controlling worker initialization profile. One of
  `"auto"`, `"local"`, or `"hpc"`.

  - `"auto"`: detects common scheduler environments (e.g. SLURM/PBS/LSF)
    and chooses `"hpc"` when detected.

  - `"local"`: standard local machine setup with default library paths.

  - `"hpc"`: cluster-oriented setup; can prepend worker library paths
    via `hpc_libpaths`.

- hpc_libpaths:

  Optional character vector of library paths to prepend on worker nodes
  when `cluster_profile = "hpc"`. Use this when compute nodes need
  explicit [`.libPaths()`](https://rdrr.io/r/base/libPaths.html) (e.g.
  non-shared user libraries). Ignored for `"local"` profile.

## Value

a S4 object of the class metime_analyzer with analysis results appended
to the result section.

## Details

Add details here
