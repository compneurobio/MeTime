# Case Study: Multi-Cohort Meta Analysis

> One-line summary: combine multiple datasets and compare consistency of findings.

## Objective

`TODO: explain cross-cohort question.`

## Cohort harmonization checklist

- Feature name normalization
- Timepoint alignment
- Covariate alignment

## Pipeline template

```r
# obj <- mod_merge_metime_analysers(...)
# obj <- calc_<method>(obj, which_data = c("cohort_a", "cohort_b"), ...)
# obj <- meta_<method>(obj, ...)
```

## Interpretation

- Concordant findings across cohorts
- Heterogeneity indicators
