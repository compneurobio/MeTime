# Calculation of generalized additive mixed models (GAMMs)

Fits multiple generalized additive mixed models (GAMMs) to a
longitudinal dataset.

## Usage

``` r
calc_gamm(
  object,
  which_data,
  stratifications = NULL,
  cols_for_meta = NULL,
  threshold = c("none", "nominal", "li", "fdr", "bonferroni"),
  verbose = T,
  name = "regression_gamm_1",
  num_cores = NULL,
  k = 10,
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

  a character of length 1 to define the type of threshold for
  significant interactions. allowed inputs are "li", "FDR", "bonferroni"
  and "nominal"(cutoff p=0.05, set as Default)

- verbose:

  a logical on whether to print the calculation progress. Default set to
  FALSE.

- name:

  a character vector to define the index within the results. Should be
  equal to length of which_data. Default set to regression_gamm_1.

- num_cores:

  numeric input to define the number of cores that you want to use for
  parallel computing. Default is set to NULL which is
  parallel::detectCores() -1.

- k:

  numeric input for setting the basis complexity for smoothing in gam.
  See more on [mgcv::bam](https://rdrr.io/pkg/mgcv/man/bam.html).
  Default is set to 10 as suggested by mgcv

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

- random:

  a character vector defining which variables should be treated as
  random effects. Default set to "subject".

- interaction:

  a character vector defining which interaction terms should be added to
  the model. Default set to NULL, with no interaction added.

## Value

a S4 object of the class metime_analyzer with analysis results appended
to the result section.

## Details

The calculation function fits multiple generalized additive mixed models
(GAMMs) on a longitudinal dataset. Here, one model fits one metabolite
vs one trait. The degree of smoothness of a model term is estimated as
part of the fitting.
