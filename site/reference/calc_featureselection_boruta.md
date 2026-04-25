# Function to calculate dependent variables

Calculate feature importance of predictors onto a response vector using
the Boruta algorithm.

## Usage

``` r
calc_featureselection_boruta(
  object,
  which_x,
  which_y,
  verbose = TRUE,
  name = "boruta_1",
  cols_for_meta_x = NULL,
  cols_for_meta_y = NULL,
  maxRuns = 100,
  num_cores = NULL,
  save_per_run = F,
  run_index = NULL,
  cluster_profile = c("auto", "local", "hpc"),
  hpc_libpaths = NULL
)
```

## Arguments

- object:

  a S4 object of the class "metime_analyzer".

- which_x:

  a character defining the name of the dataset containing the
  predictors.

- which_y:

  a character defining the name of the dataset containing the response.
  Can be factor for classification or numeric vector for regression.

- verbose:

  a logical whether the progress should be printed. Default set to TRUE.

- cols_for_meta_x:

  A named list of a character vector to define columns for meta info of
  which_x Names of the character vector in this case should always start
  with "color\_".

- cols_for_meta_y:

  a named list of a character vector to define columns for meta info of
  which_y

- num_cores:

  numeric input to define the number of cores that you want to use for
  parallel computing. Default is set to NULL which is
  parallel::detectCores() -1.

- save_per_run:

  (Experimental) a logical on whether results should be saved in the
  same directory/tmp as single csv.

- run_index:

  (Experimental) a vector of runs corresponding the row number of a

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

List of conservation index results

## Details

The Boruta algorithm is a feature selection method that compares the
importance of each variable to that of randomly generated shadow
variables. It uses an iterative two step approach:

1.  the algorithm creates shadow variables that are random permutations
    of the original variables.

2.  Secondly, it fits a machine learning model to the original and
    shadow variables, and calculates the importance of each variable
    based on the difference between the model performance with and
    without the variable.

Variables significantly more important than their shadow counterparts
are considered important and selected for the final model. This process
is repeated until all variables have been evaluated. The Boruta
algorithm is particularly useful for datasets with many variables and
noisy data, and can help improve the accuracy and interpretability of
predictive models.

## See also

[Boruta::Boruta](https://rdrr.io/pkg/Boruta/man/Boruta.html)
