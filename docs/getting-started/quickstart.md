# Quickstart

> Run your first MeTime pipeline with bundled sample data.

## Goal

In this guide, you will run a minimal pipeline for dimensionality reduction.

## Minimal runnable example with HuMet data

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

## Next

- Continue to [`data_preparation.md`](data_preparation.md).
