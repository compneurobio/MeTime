# Testing and Validation

> One-line summary: minimum checks before merging new methods or docs.

## Unit test checklist

- Input validation paths
- Expected output structure
- Error/warning behavior
- Reproducibility (`set.seed`)

## Documentation test checklist

- Every code chunk runs
- Every link resolves
- Every page has troubleshooting notes

## Recommended CI checks

```bash
R CMD check .
Rscript -e "pkgdown::build_site(preview = FALSE)"
```
