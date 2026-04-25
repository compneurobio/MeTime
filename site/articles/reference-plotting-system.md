# Plotting System

## Plotting System

> Generate and store plots from MeTime results.

### Two plotting pathways

1.  `plot(<metime_analyser>)` for direct rendering.
2.  [`mod_generate_plots()`](https://compneurobio.github.io/MeTime/reference/mod_generate_plots.md)
    to generate and optionally store plot objects in `results$plots`.

Supported `type` values include:

- `network`
- `PCA`
- `UMAP`
- `tSNE`
- `CI_metabotype`
- `CI_metabolite`
- `regression`
- `pairwise_distance`
- `pairwise_correlation`
- `colinearity`
- `distribution_samples`
- `distribution_metabs`
- `feature_selection`
- `clusters`

### Plot template

``` r
# p <- plot(obj, results_index = "<result>", ...)
# obj <- add_plot(obj, name = "<plot_name>", plot = p)
# obj <- mod_generate_plots(obj, results_index = "<result>", .interactive = FALSE, type = "<type>")
```

### Plot QA checklist

- Axis labels are interpretable.
- Group/color encodings match metadata.
- Legends and units are present.
