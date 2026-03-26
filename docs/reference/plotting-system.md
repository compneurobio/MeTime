# Plotting System

> One-line summary: generate and store plots from MeTime results.

## Two plotting pathways

1. `plot(<metime_analyser>)` for direct rendering.
2. `mod_generate_plots()` for storing all possible plots for a paricular result in object.
In mod_generate_plots, <type> argument is used to select the result you want to analyse the follwing inputs are allowed: "network", "PCA", "UMAP", "tSNE", "CI_metabotype", "CI_metabolite", "regression", "pairwise_distance", "pairwise_correlation", "colinearity", "regression", "distribution_samples, "distribution_metabs", "feature_selection", "clusters" 

## Plot template

```r
# p <- plot(obj, results_index = "<result>", ...)
# obj <- add_plot(obj, name = "<plot_name>", plot = p)
# obj <- mod_generate_plots(obj, results_index="<result>", .interactive=FALSE, type="<type>")
```

## Plot QA checklist

- Axis labels are interpretable.
- Group/color encodings match metadata.
- Legend and units are present.
