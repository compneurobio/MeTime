# Data preparation 
> One-line summary: Create your own analyser object for analysing

## Goal

In this guide, you will:

1. Build a `metime_analyser` object.
2. Apply one `mod_*` transformation.
3. Run one `calc_*` analysis.
4. Create one plot.

## Building your own analyser object

This id a general package that can handle any type of longitudinal dataset and we expect users to make the following changes to the dataset
1. The sample ids should always be in this format: [a-z|A-z][0-9]+_[a-z|A-Z][0-9]+ (Example: subject=R1, time=t0, id=R1_t0). The part before the underscore represents the subject and the part after represents the timepoint of measurement. If the timepoints in the data are not a singular value then we suggest the user to create a psuedo timescale to match this format.
2. Every row_data dataframe should contain the columns id, subject and time and every col_data dataframe should contain the column id. And the ids in row_data should match the rownames of the data matrix and the ids in col_data should match the colnames of the data matrix.

Once the data has been prepared in the following format, you can build your own analyser object in the following ways

### Loading all the data at once from a folder with all files

The user can store all the files that contain the data in a directory and can parse the path to create an object. However, it is important to note that the naming of the files(extension should either be .rds or .RDS) for a particular dataset "test" should be in this pattern:
 data: test_data
 col_data: test_col_data
arow_data: test_row_data

```r
path <- "/path/to/directory"

annotations=list(phenotype="name_of_the_dataset", 
                medication="name_of_the_dataset") # this applies only if you have multiple omics dataset from same individuals else set to NULL
object <- get_files_and_names(path, annotations_index=annotations)

```

Or you can add datasets manually using the following format:

```r

data <- read.csv("your/data/file.csv")
row_data <- read.csv("your/row_data/file.csv")
col_data <- read.csv("yout/col_data/file.csv")

object <- get_make_analyser_object(data=data.frame, # dataframe containing data
                                  col_data=data.frame, # dataframe containing col_data
                                  row_data=data.frame, # dataframe containing row_data
                                  annotations_index=annotations,
                                  name="name of the dataset")

# One can also add datasets afterwards in this way:
object <- object %>% add_dataset(data=data, row_data=row_data, col_data=col_data, name="name")

```

## If this fails for your data

- run validate_metime_analyser()
- Check sample IDs format.
- Check if required metadata columns exist such as id, subject and/or time.

## Applying mod and calc functions - example to perform conservation index analysis

```r
# update which_data variable as per your dataset name
which_data <- "name"
# Using mod_mutate (wrapper of dplyr::mutate for ease) and running PCA
object <- object %>% 
            mod_mutate(type="row_data", which_data=which_data,
                    test_col=seq_along(rownames(get_rowdata(object, which_data=which_data)))) %>%
            mod_trans_zscore(which_data=which_data) %>%
            calc_conservation_metabotype(which_data=which_data, 
                               stratifications=list(time=c("t0", "t12", "t24")),
                               verbose=F, 
                               cols_for_meta=NULL, 
                               name="Conservation_Metabotype_scaled")
# if you want plots to be colored by other variables add cols_for_meta as a character vector with colnames from row_data in the above example

# generate plots
object <- object %>% mod_generate_plots(type="CI_metabotype", .interactive=T)

# generate report

object %>% write_report(file="dimensionality_reduction.html", title="PCA results")

```


## Next

- Continue to [`sample-data-walkthrough.md`](sample-data-walkthrough.md).
