# MeTime

MeTime (**Metabolomics Time**) is an R package for modular, reproducible analysis of longitudinal metabolomics data.

## What MeTime provides

- A central S4 object (`metime_analyser`) to manage multiple datasets and metadata.
- Pipeline-style functions:
  - `mod_*` for data preparation/transformation
  - `calc_*` for analyses
  - `meta_*` for cross-result aggregation
- Built-in support for common longitudinal workflows:
  - distributions and distances
  - dimensionality reduction
  - imputation and feature selection
  - conservation index
  - regression models
  - network analysis

## Install

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}

options(repos = BiocManager::repositories())
devtools::install_github("compneurobio/MeTime", dependencies = TRUE)
```

## Quick example

```r
library(MeTime)
data("humet_object")

humet_object <- humet_object %>%
  mod_trans_zscore(which_data = "humet_subset_data") %>%
  calc_dimensionality_reduction_samples(
    which_data = "humet_subset_data",
    type = "PCA",
    cols_for_meta = c("Factor.Challenge.Value.", "Factor.Challenge.Value.Day."),
    name = "PCA_samples"
  ) %>%
  mod_generate_plots(type = "PCA", .interactive = TRUE)
```

## Documentation

- Docs hub: [`docs/README.md`](docs/README.md)
- Getting started: [`docs/getting-started/`](docs/getting-started)
- User guides: [`docs/user-guides/`](docs/user-guides)
- Reference guides: [`docs/reference/`](docs/reference)
- Case studies: [`docs/case-studies/`](docs/case-studies)
- Package reference (`?function_name`) and `vignettes/`

## Website

A pkgdown site is built from repository docs/vignettes via CI using `.github/workflows/pkgdown.yaml`.
