# Quickstart (10 minutes)

> One-line summary: run your first MeTime pipeline with sample data.

## Goal

In this guide, you will:

1. Build a `metime_analyser` object.
2. Apply one `mod_*` transformation.
3. Run one `calc_*` analysis.
4. Create one plot.

## Minimal runnable example

```r
library(MeTime)

# TODO: replace with real package dataset and function names used in the package
# data(humet_data)
# obj <- get_make_analyser_object(...)
# obj <- mod_<name>(obj, ...)
# obj <- calc_<name>(obj, ...)
# p <- plot(obj, results_index = "<result_name>")
# p
```

## Expected result

- Output object has at least one named entry in `@results`.
- A plot is rendered without errors.

## If this failed

- Check sample IDs format.
- Check required metadata columns exist.
- Validate dataset names in `which_data`.

## Next

- Continue to [`sample-data-walkthrough.md`](sample-data-walkthrough.md).
