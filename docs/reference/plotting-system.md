# Plotting System

> One-line summary: generate and store plots from MeTime results.

## Two plotting pathways

1. `plot(<metime_analyser>)` for direct rendering.
2. `mod_generate_plots()` for storing all possible plots for a paricular result in object.

## Plot template

```r
# p <- plot(obj, results_index = "<result>", ...)
# obj <- add_plot(obj, name = "<plot_name>", plot = p)
```

## Plot QA checklist

- Axis labels are interpretable.
- Group/color encodings match metadata.
- Legend and units are present.
