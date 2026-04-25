# Troubleshooting

## Troubleshooting

> One-line summary: common errors, root causes, and fixes.

This page includes practical checks you can run quickly with the bundled
`humet_object` example dataset. Some sections include placeholders that
you can later replace with project-specific details.

### Quick diagnostics first

Run this minimal diagnostic block before deep debugging:

``` r
library(MeTime)
data("humet_object")

# 1) Object validity
validate_metime_analyser(humet_object)

# 2) Inspect available datasets
names(humet_object@list_of_data)

# 3) Inspect available results
names(humet_object@results)

# 4) Basic summary
summarize_dataset(humet_object, which_data = "humet_subset_data")
```

### Top issues

#### 1) Dataset not found in object

- **Symptom:** `which_data` errors or empty operations.
- **Likely cause:** typo or missing dataset name.
- **Quick fix:**

``` r
names(humet_object@list_of_data)
# then use one of these names in which_data
```

#### 2) Plot generation fails

- **Symptom:** no plot output or plotting error.
- **Likely cause:** missing required result type / metadata fields.
- **Quick fix:**

``` r
# inspect result metadata
humet_object@results[[1]]$information$calc_type

# then pick matching type in mod_generate_plots(...)
```

#### 3) Invalid sample IDs

- **Symptom:** merge/filter steps produce unexpected row counts.

- **Likely cause:** inconsistent sample ID format across datasets.

- **Quick fix:** standardize IDs before building/merging data.

- Canonical ID regex: `[a-z][A-Z][0-9]+_[a-z][A-Z][0-9]+`

#### 4) Too many missing values (NAs)

- **Symptom:** analyses return empty or unstable results.
- **Likely cause:** high missingness at feature or sample level.
- **Quick fix:**

``` r
humet_object %>%
  mod_filter_features_by_missingness(which_data = "humet_subset_data", threshold = 0.3) %>%
  mod_filter_samples_by_missingness(which_data = "humet_subset_data", threshold = 0.3)
```

#### 5) Model formula / covariate mismatch

- **Symptom:** regression functions fail or drop terms.
- **Likely cause:** covariates not present, wrong type (numeric/factor),
  or invalid formula columns.
- **Quick fix:**

``` r
# inspect columns before running model
colnames(humet_object@list_of_data[["humet_subset_data"]])
```

#### 6) Package/dependency install issues

- **Symptom:** install fails due to system libraries.
- **Likely cause:** missing OS dependencies for compiled packages.
- **Quick fix:** install listed system libraries from
  `docs/getting-started/installation.md`.

### Recommended debug workflow

1.  Reproduce with `humet_object` first to isolate whether issue is
    data-specific.
2.  Run pipeline one step at a time and inspect intermediate object
    states.
3.  Use
    [`validate_metime_analyser()`](https://compneurobio.github.io/MeTime/reference/validate_metime_analyser.md)
    after each major `mod_*` step.
4.  Save intermediate objects with
    [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) to enable
    reproducible bug reports.
5.  Confirm that transformations match intended scale/log settings.

### See fixed bugs

Check closed GitHub issues for solved examples and known edge cases.
