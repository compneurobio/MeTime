# Data preparation

> Create a `metime_analyser` object from your own data.

## Goal

In this guide, you will:

1. Build a `metime_analyser` object.
2. Apply one `mod_*` transformation.
3. Run one `calc_*` analysis.
4. Create one plot.

## Required data conventions

1. Sample IDs should follow `<subject>_<timepoint>` (example: `R1_t0`) and should be readable by this regex format: "[a-z][A-z][0-9]+_[a-z][A-Z][0-9]+".
2. Every `row_data` table must include `id`, `subject`, and `time`.
3. Every `col_data` table must include `id`.
4. `rownames(data)` must match `row_data$id`, and `colnames(data)` must match `col_data$id`.
5. All tables should be of data.frame() class

## Build an analyser object

### Load all files from a directory

For a dataset called `test`, expected file names are:
- `test_data.rds`
- `test_col_data.rds`
- `test_row_data.rds`

Annotations are currently built for cohort data that is multi-platform. For example if a longitudinal cohort is measured on different metabolomics platforms this would lead to multiple datasets, however the clinical information or the medication intake information remains the same for all. In order to remove duplications (adding this to each dataset's row_data) we provide a separate holder for such data and we expect you to annotate them. 
See regression analysis in user-guides/building-pipelines.md for more information on how to handle such data

In this extended example for a cohort A, expected file names are:
- Data files:
  - `platform1_data.rds`
  - `platform2_data.rds`
  - `phenotype_data.rds`
  - `medication_data.rds`
- And the respective row_data and col_data files as described above

```r
path <- "/path/to/directory"

annotations <- list(
  phenotype = "phenotype_data", # can be set to NULL if not applicable
  medication = "medication_data" # can be set to NULL if not applicable
) 
object <- get_files_and_names(path, annotations_index = annotations)
```

### Create object from in-memory data frames

```r
data <- read.csv("your/data/file.csv")
row_data <- read.csv("your/row_data/file.csv")
col_data <- read.csv("your/col_data/file.csv")

object <- get_make_analyser_object(
  data = data,
  col_data = col_data,
  row_data = row_data,
  annotations_index = annotations,
  name = "name_of_the_dataset"
)

# Add another dataset later
object <- object %>%
  add_dataset(data = data, row_data = row_data, col_data = col_data, name = "name_2")
```

## Example: conservation index pipeline

```r
which_data <- "name"

object <- object %>%
  mod_mutate(
    type = "row_data",
    which_data = which_data,
    test_col = seq_along(rownames(get_rowdata(object, which_data = which_data)))
  ) %>%
  mod_trans_zscore(which_data = which_data) %>%
  calc_conservation_metabotype(
    which_data = which_data,
    stratifications = list(time = c("t0", "t12", "t24")), # change accordingly
    verbose = FALSE,
    cols_for_meta = NULL,
    name = "Conservation_Metabotype_scaled"
  ) %>%
  mod_generate_plots(type = "CI_metabotype", .interactive = TRUE)

object %>% write_report(file = "CI_results.html", title = "Conservation index analysis")
```

## Merge analyser objects

```r
annotations <- object1@annotations
object <- mod_merge_metime_analysers(
  object1,
  object2,
  annotations_index = annotations
)
```

## Next

- Continue to the user guides for analysis-specific pipelines.
