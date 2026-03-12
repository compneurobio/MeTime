# Troubleshooting

> One-line summary: common errors, root causes, and fixes.

## Top issues

### 1) Dataset not found in object
- **Cause:** `which_data` does not match dataset names.
- **Fix:** check `names(object@list_of_data)`.

### 2) Plot generation fails
- **Cause:** incompatible result type or missing metadata.
- **Fix:** inspect `object@results[[<index>]]$information$calc_type`.

### 3) Invalid sample IDs
- **Cause:** IDs not in required format.
- **Fix:** standardize IDs before object creation.

### 4) See fixed bugs here

## Debug workflow template

```r
# print names and structures before each calc step
```
