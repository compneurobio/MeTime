# An automated fucntion to calculate GGM from multibipartite lasso approach

automated funtion that can be applied on s4 object of class
metime_analyser to calculate a network using multibipartite lasso
Conceptualization & Methodology: Richa Batra and Jan Krumsiek.
Implementation: Richa Batra, Bharadwaj Marella, Patrick Weinisch and
Matthias Arnold.

## Usage

``` r
calc_ggm_multibipartite_lasso(
  object,
  which_data,
  alpha = 1,
  nfolds = 3,
  stratifications,
  cols_for_meta,
  cluster_profile = c("auto", "local", "hpc"),
  num_cores = NULL,
  hpc_libpaths = NULL,
  name = "mbpl_ggm_1",
  ...
)
```

## Arguments

- object:

  S4 object of cĺass metime_analyser

- which_data:

  a character or a character vector naming the datasets of interest

- alpha:

  tuning parameter for lasso + ridge regression in glmnet. Default set
  to 1 to perform LASSO

- nfolds:

  nfolds for cv.glmnet. Default set to 3

- stratifications:

  List to stratify data into a subset. Usage list(name=value)

- cols_for_meta:

  a list of character vectors of column names to be used for
  visualization of the networks.

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

- name:

  character to define the name of the results

- ...:

  additional arguments for cv.glmnet function

## Value

Analyser object with updated results of this calculation
