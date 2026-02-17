# Building Pipelines

> One-line summary: chain modular functions into reproducible analysis pipelines.

## Pipeline mental model

1. Build object
2. Modify data (`mod_*`)
3. Calculate (`calc_*`)
4. Meta-analysis (`meta_*`, optional)
5. Plot/export

## Starter pipeline template

```r
library(MeTime)

# obj <- get_make_analyser_object(...)
# obj <- mod_<step1>(obj, ...)
# obj <- mod_<step2>(obj, ...)
# obj <- calc_<analysis>(obj, which_data = "<dataset>", name = "<result_1>")
# obj <- add_plot(obj, ...)
# write_results(obj, ...)
```

## Reproducibility checklist

- Save package version
- Store parameters and seed
- Export result tables and plots
