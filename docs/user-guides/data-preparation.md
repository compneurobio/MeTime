# Data Preparation

> One-line summary: make longitudinal datasets compatible with MeTime.

## Required input shape

- Samples in rows and metabolites/features in columns (or specify package convention).
- Consistent dataset naming across object slots.

## Required ID format

- Pattern: `<regex>`
- Example valid IDs: `<examples>`
- Example invalid IDs: `<examples>`

## Required metadata columns

| Column | Description | Required | Example |
|---|---|---:|---|
| `<col>` | `<description>` | yes | `<value>` |

## Validation snippet

```r
# TODO: add validation helper calls
```

## Common pitfalls

- Mismatched sample IDs across matrices.
- Missing timepoint encoding.
- Non-numeric concentration fields.
