# Testing and Validation

## Testing and Validation

> One-line summary: minimum checks before merging new methods or docs.

This page provides a lightweight testing framework with placeholders you
can refine later.

### Test strategy

Use a layered approach:

1.  **Smoke tests**: package loads and example dataset runs.
2.  **Unit tests**: input validation + output contract checks.
3.  **Integration tests**: short end-to-end pipeline on `humet_object`.
4.  **Documentation checks**: examples and links remain valid.

### Quick local checks

``` bash
R CMD build .
R CMD check --no-manual --as-cran MeTime_*.tar.gz
Rscript -e "pkgdown::build_site(preview = FALSE)"
```

### Suggested testthat structure (placeholder)

``` text
tests/
  testthat/
    test-smoke.R
    test-modules.R
    test-calculations.R
    test-regressions.R
    test-plotting.R
```

- `<PLACEHOLDER_ADDITIONAL_TEST_FILES>`
- `<PLACEHOLDER_PERFORMANCE_TEST_PLAN>`

### Minimum test checklist

Input validation paths covered for each new public function.

Expected output structure verified (`class`, key columns, metadata).

Error/warning behavior verified (`expect_error`, `expect_warning`).

Reproducibility verified with
[`set.seed()`](https://rdrr.io/r/base/Random.html) where relevant.

New docs examples run without manual intervention.

### Example smoke test (humet_object)

``` r

library(testthat)
library(MeTime)

test_that("example dataset pipeline runs", {
  data("humet_object")

  out <- humet_object %>%
    mod_trans_zscore(which_data = "humet_subset_data") %>%
    calc_dimensionality_reduction_samples(
      which_data = "humet_subset_data",
      type = "PCA",
      cols_for_meta = c("Factor.Challenge.Value.", "Factor.Challenge.Value.Day."),
      name = "PCA_samples"
    )

  expect_s4_class(out, "metime_analyser")
  expect_true(length(out@results) >= 1)
})
```

### CI recommendation

Use two pipelines:

- **PR CI:** fast checks (smoke + selected unit tests).
- **Release CI:** full `R CMD check`, docs build, and artifact
  publishing.

#### Placeholder CI policy

- Required checks for merge: `<PLACEHOLDER_REQUIRED_CHECKS>`
- Allowed flaky checks: `<PLACEHOLDER_FLAKY_CHECK_POLICY>`
- Maximum CI runtime target: `<PLACEHOLDER_MAX_RUNTIME>`

### Coverage goals (placeholder)

- Global line coverage target: `<PLACEHOLDER_LINE_COVERAGE_TARGET>`
- Critical module coverage target:
  `<PLACEHOLDER_CRITICAL_PATH_COVERAGE>`
- Regression test SLA for bugs: `<PLACEHOLDER_REGRESSION_TEST_SLA>`
