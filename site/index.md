# MeTime

MeTime (**Metabolomics Time**) is an R package for modular, reproducible
analysis of longitudinal metabolomics data.

## What MeTime provides

- A central S4 object (`metime_analyser`) to manage multiple datasets
  and metadata.
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

``` r

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

``` r

library(MeTime)
data("humet_object")

humet_object <- humet_object %>%
  mod_trans_zscore(which_data = "humet_subset_data") %>%
  mod_mutate(which_data="humet_subset_data", type="row_data", 
      Factor.Value.Challenge.Day.=factor(Factor.Value.Challenge.Day.)) %>%
  calc_dimensionality_reduction_samples(
    which_data = "humet_subset_data",
    type = "PCA",
    cols_for_meta = c("Factor.Value.Challenge.", "Factor.Value.Challenge.Day."),
    name = "PCA_samples"
  ) %>%
  mod_generate_plots(type = "PCA", .interactive = TRUE)
```

## Documentation

- Docs hub:
  [`docs/README.md`](https://compneurobio.github.io/MeTime/docs/README.md)
- Getting started:
  [`docs/getting-started/README.md`](https://compneurobio.github.io/MeTime/docs/getting-started/README.md)
- User guides:
  [`docs/user-guides/README.md`](https://compneurobio.github.io/MeTime/docs/user-guides/README.md)
- Reference guides:
  [`docs/reference/README.md`](https://compneurobio.github.io/MeTime/docs/reference/README.md)
- Case studies:
  [`docs/case-studies/README.md`](https://compneurobio.github.io/MeTime/docs/case-studies/README.md)
- Developer docs:
  [`docs/developer/README.md`](https://compneurobio.github.io/MeTime/docs/developer/README.md)
- Package reference (`?function_name`)

## Citation

If MeTime contributes to your analysis or publication, please cite the
package: *Marella, B., Weinisch, P.,* Vehovec, L., Tran V., Bless, J.
J., Njipouombe Nsangou, Y. A., Kastenmüller G. & Arnold M. (2026).
MeTime: An R package for reproducible longitudinal metabolomics data
analysis. aRxiv.

For more details and the publication please see:
[`CITATION.md`](https://compneurobio.github.io/MeTime/CITATION.md)

## Website

A pkgdown site is built from repository docs/vignettes via CI using
`.github/workflows/pkgdown.yaml`.

## Container image

A Docker image can be published to GitHub Container Registry (GHCR) via
`.github/workflows/docker-publish.yaml` when pushing release tags (for
example `v1.0.1`).
