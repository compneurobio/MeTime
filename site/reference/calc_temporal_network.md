# An automated function to caluclate temporal network with lagged model

calculates temporal networks for each dataset with a lagged model as
used in graphical VAR

## Usage

``` r
calc_temporal_network(
  object,
  which_data,
  lag,
  stratifications,
  alpha = 1,
  nfolds = 3,
  cols_for_meta,
  names,
  cluster_profile = c("auto", "local", "hpc"),
  num_cores = NULL,
  hpc_libpaths = NULL,
  ...
)
```

## Arguments

- object:

  S4 object of class metime_analyser

- which_data:

  dataset or datasets to be used

- lag:

  which lagged model to use. 1 means one-lagged model, similarly
  2,3,..etc

- stratifications:

  List to stratify data into a subset. Usage list(name=value)

- alpha:

  parameter for regression coefficient. Set to 1 for lasso regression

- nfolds:

  nfolds parameter for glmnet style of regression. Default is set to 3

- cols_for_meta:

  a list of character vectors of column names to be used for
  visualization of the networks.

- names:

  character vector with the same length as that of possible models

- cluster_profile:

  Character string controlling worker initialization profile. One of
  `"auto"`, `"local"`, or `"hpc"`.

  - `"auto"`: detects common scheduler environments (e.g. SLURM/PBS/LSF)
    and chooses `"hpc"` when detected.

  - `"local"`: standard local machine setup with default library paths.

  - `"hpc"`: cluster-oriented setup; can prepend worker library paths
    via `hpc_libpaths`.

- num_cores:

  numeric input to define the number of cores that you want to use for
  parallel computing. Default is set to NULL which is
  parallel::detectCores() -1.

- hpc_libpaths:

  Optional character vector of library paths to prepend on worker nodes
  when `cluster_profile = "hpc"`. Use this when compute nodes need
  explicit [`.libPaths()`](https://rdrr.io/r/base/libPaths.html) (e.g.
  non-shared user libraries). Ignored for `"local"` profile.

- ...:

  additional arguments for cv.glmnet function

## Value

S4 object with updated temporal network results
