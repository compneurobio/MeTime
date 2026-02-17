# metime_analyser Class

> One-line summary: understand object slots, invariants, and lifecycle.

## Slots

| Slot | Type | Description |
|---|---|---|
| `list_of_data` | list | concentration matrices |
| `list_of_row_data` | list | sample metadata |
| `list_of_col_data` | list | feature metadata |
| `annotations` | list | phenotype/medication definitions |
| `results` | list | functions, data, plots, information |

## Invariants

- Dataset names must align across list slots.
- Result names should be unique.
- Metadata should remain synchronized after filtering.

## Inspection helpers

```r
# show(obj)
# get_data(obj, which_data = "...")
# get_results(obj)
```
