# Troubleshooting

> One-line summary: common errors, root causes, and fixes.

## Top issues

### 1) Dataset not found in object
- **Cause:** `which_data` does not match dataset names.
- **Fix:** check `names(object@list_of_data)`.

### 2) Plot generation fails
- **Cause:** incompatible result type or missing metadata.
- **Fix:** inspect `object@results[[<index>]]$information$calc_type` and choose the relevant result_type or plot_type.

### 3) Invalid sample IDs
- **Cause:** IDs not in required format.
- **Fix:** standardize IDs before object creation.

### 4) See fixed bugs

- check the closed issue page to see if a similar problem has been previously solved

## Debug workflow ideas:

1. Run the pipeline step by step and run validate_metime_analyser() to check where the issue lies
2. In case everything works and the result is empty - check if data is has no NAs
3. Use Codex or similar AI tool to inspect


