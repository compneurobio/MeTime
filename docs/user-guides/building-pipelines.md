# Building Pipelines

> One-line summary: chain modular functions into reproducible analysis pipelines.

## Pipeline mental model

1. Build object
2. Modify data (`mod_*`)
3. Calculate (`calc_*`)
4. Plot/export 
4. Meta-analysis (`meta_*`, optional)

## Starter pipeline template ()

```r
library(MeTime)

# obj <- get_make_analyser_object(...)
# obj <- mod_<step1>(...) %>%
#        mod_<step2>(...) %>%
#        calc_<analysis>(which_data = "<dataset>", name = "<result_1>") %>% 
#        mod_generate_plots(...) %>%
#        write_report(...)
```

## Reproducibility checklist

- Save package version
- Store parameters and seed
- Export result tables and plots
