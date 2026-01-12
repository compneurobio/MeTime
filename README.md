# Introduction to MeTime

MeTime (Metabolomics Time) helps scientists analyze longitudinal metabolomics data without rebuilding bespoke pipelines for every study. The package provides modular, pipeline-friendly functions that allow you to chain analyses for a given dataset, then reuse or remix those steps across projects. Example workflows are provided as vignettes (see the **Vignettes** section below). 

To make this package user-friendly and easy to implement, we use a puzzle analogy. Each analysis answers a specific question and can be considered as a sub-puzzle. These sub-puzzles collectively give a complete mechanistic picture, and each sub-puzzle is made up of modular pipeline functions. See the images below for a visual depiction of this approach. To understand the modular pieces shown in the image, please refer to section [2. Building pipelines].

<img width="463" alt="puzzle" src="https://user-images.githubusercontent.com/64539275/232745515-a4bfc9fe-d353-402a-92b5-7433b0aed2a4.PNG">

There are different methods to analyse longitudinal metabolomics data such that each method has its own unique significance and answers a specific question. These methods are listed below. To better understand the application of each of these methods and to see examples on how to perform such an analysis, please refer to their specific vignettes.
<br> 1.  Distributions  
<br> 2.  Feature selection
<br> 3.  Imputation
<br> 4.  Dimensionality reduction
<br> 5.  Eigendata calculation
<br> 6.  Conservation index analysis
<br> 7.  Regressions
<br> 8.  Data-driven networks
<br> 9.  meta-analyses


## Vignettes

Each vignette mirrors a master script and focuses on a single analysis workflow:

1. Basic data statistics  
2. Feature selection  
3. Imputation  
4. Dimensionality reduction  
5. Eigendata calculation  
6. Conservation index analysis  
7. Regressions  
8. Data-driven networks  

## Getting started

### Installation

```r
install.packages("remotes")
remotes::install_github("your-org/MeTime")
```

### Quick start

```r
library(MeTime)
library(magrittr)

# Load or build a metime_analyser object
# object <- get_make_analyser_object(...)

# Example pipeline
object %>%
  mod_trans_zscore(which_data = "nmr_data") %>%
  calc_dimensionality_reduction_samples(which_data = "nmr_data", type = "PCA")
```

### Building the website (pkgdown)

```r
pkgload::load_all()
pkgdown::build_site()
```

## Documentation map

This document is divided into four main sections: 
<br> [1. metime_analyser class and data preparation]
<br> [2. Building pipelines]
<br> [3. For Developers]

## 1. metime_analyser class and data preparation

### 1.1. metime_analyser class

This package builds upon the S4 class of metime_analyser which serves as a central object that stores the data, results and the information of functions applied onto this object in a pipeline. The reasons to create such an object is as follows:
<br> 1. There are often multiple datasets that a user wants to analyse and there is no class in R that can store all the data at once. The closest relative of metime_analyser is summarizedExperiment(SE), which however stores information of only a single dataset that is cross-sectional in nature
<br> 2. As the metime_analyser object will contain multiple datasets at once. It is easier to parse this into other functions and modify/analyse all the datasets at once thereby removing the need of duplicating the same analysis for different datasets.
<br> 3. Moreover, users can now perform analyses which clubs two or more datasets at once.
<br> 4. To be able to reproduce the results and maintain transparency, information regarding the functions applied onto a metime_analyser object are also stored. See the structure of results below to understand this better. 

The metime_analyser class has 5 slots:
<br> 1. list_of_data: Consists of a list of data matrices of metabolite concentrations
<br> 2. list_of_row_data: Consists of a list of row-data information(samples) for the respective data matrices
<br> 3. list_of_col_data: Consists of a list of col-data information(metabolites) for the respective data matrices
<br> 4. annotations: This is a list to define how the phenotype data and medication data are named. Can also include other datasets, however, only phenotype and medication data are important as they are different from the other datasets which are actually analysed.
<br> 5. results: List where the results of the analyses are stored. This list contains up to 4 elements namely:
<br> &emsp; 1. functions_applied: A named list with names being the functions that were applied onto the object until a particular calculation(analysis) is performed and the values of each list item is a named character vector with parameters as names of the character vector and values are their respective values. If an argument's value is anything other than a character vector then that value is converted into a class of character type. 
 <br> &emsp; 2. plot_data: A list of dataframes which are the results of an analysis.
 <br> &emsp; 3. plots: A list where the plots generated by mod_generate_plots() are stored. For more information see [Building pipelines] section
 <br> &emsp; 4. information: A list with two elements calc_type and calc_info. calc_type is a character vector that defines the type of analysis and calc_info is a more detailed description of the calculation performed. calc_type and calc_info are generated automatically after a calculation and is used as key for plotting the results or for performing meta analysis. 

<img width="435" alt="metime_analyser" src="https://user-images.githubusercontent.com/64539275/232745616-7149a04a-a051-4fea-9079-5f5fb2106089.png">

<img width="376" alt="metime_analyser_results" src="https://user-images.githubusercontent.com/64539275/232745659-0c9fcea9-145f-4829-881a-f279d035d79f.png">

### 1.2. Data preparation

This package aims to be a general package that can handle any type of longitudinal dataset and in order to achieve that we expect the users to make a few changes to the dataset. These changes are:
1. The sample ids should always be in this format: [a-z|A-z][0-9]+_[a-z|A-Z][0-9]+(Example: subject=R1, time=t0, id=R1_t0). The part before the underscore represents the subject and the part after represents the timepoint of measurement. If the timepoints in the data are not a singular value then we suggest the user to create a pseudotime scale to match this format.
2. Every row_data dataframe should contain the columns id, subject and time and every col_data dataframe should contain the column id. And the ids in row_data should match the rownames of the data matrix and the ids in col_data should match the colnames of the data matrix. 

The first step in using this package is to create an S4 object of class metime_analyser with all the data that a user wants to analyse. There are multiple ways in which this S4 object can be created:

#### 1.2.1. Loading all the data at once from a folder with all files

The user can store all the files that contain the data in a directory and can parse the path to create an object. However, it is important to note that the naming of the files(extension should either be .rds or .RDS) for a particular dataset "test" should be in this pattern:
 data: test_data
 col_data: test_col_data
arow_data: test_row_data

```{r, eval=FALSE}
path <- "/path/to/directory"
annotations=list(phenotype="name_of_the_dataset/file", 
                      medication="name_of_the_dataset/file")
object <- get_files_and_names(path, annotations_index=annotations)
```
    
#### 1.2.2 Using data-frames to make your own object 

This option is for users who want to analyse only one dataset at once. 

```{r, eval=FALSE}
annotations=list(phenotype="name_of_the_dataset/file", 
                      medication="name_of_the_dataset/file")
object <- get_make_analyser_object(data=data.frame, # dataframe containing data
                                  col_data=data.frame, # dataframe containing col_data
                                  row_data=data.frame, # dataframe containing row_data
                                  annotations_index=annotations,
                                  name="name of the dataset")
```

#### 1.2.3. Appending a new dataset to a previously created object 

```{r, eval=FALSE}
object <- get_append_analyser_object(object, data, col_data, row_data, name)
```

#### 1.2.4. merge two or more analyser objects

```{r, eval=FALSE}
annotations=object1@annotations
object <- mod_merge_metime_analysers(object1, object2, ...,
                                      annotations_index=annotations)
                                      
```

## 2. Building pipelines

To better understand how to build your own pipelines, it is important to understand the different types of functions that are available in this package.

### 2.1. Functions available in this package 

The functions of this package can be vastly categorized into two kinds: pipeline functions and non-pipeline functions. Pipeline functions are functions that are used in a pipeline and that take in metime_analyser object as an input and return a metime_analyser object as output. Non-pipeline functions are functions that either take in metime_analyser object as an input or return metime_analyser object but not both. One can imagine the non-pipeline functions to be either at the starting of a pipeline or at the end of a pipeline but never in-between. 

Moreover, each function comes with a prefix that assists the user in understanding the role of this function. The different types of prefixes and their roles are as follows:

```{r functions, eval=F}
knitr::include_graphics("vignettes/images/functions.png")
```

### 2.2. Pipelines and analysis

Pipeline functions take in object as a primary argument and returns the object with modifications based on the function. This, hence, enables users to to work with a pipe operator, either the |> operator in base R or the %>% operator from the [magrittr](https://cran.r-project.org/web/packages/magrittr/vignettes/magrittr.html) package. The pipe operator eliminates the need for repeated assignments or temporary variables to provide smooth connections between pipeline steps to create pipelines that are straightforward to use, highly modular, and reproducible. To see examples refer to vignettes of methods listed above.
  
### 2.3. Plotting with MeTime

Visualizing results the right way is extremely important to better grasp the biological insights. To this end, we have developed two functions that will enable users to quickly generate relevant plots specific to a method. One function is a pipeline function that can be used in the pipeline in order to let the user store the plots in the metime_analyser object and the other is a standalone S4 method of plot() for metime_analyser class. The former function is a wrapper of the latter function and the latter function is a wrapper of ggplot/visNetwork based on the plot type and the additional arguments in plot() function are for adding aesthetics to the plot. Below is a pictorial representation of the two methods and how they are different. 

![plotting](https://user-images.githubusercontent.com/64539275/232745864-d2d87a11-2648-497e-a2c7-80efadf30d2e.png)

## 3. For developers

As mentioned previously all MeTime functions can be classified into pipeline functions and non-pipeline functions. Here we describe the skeleton functions for each of these and the essence behind the naming. Moreover all pipeline functions are S4 methods in nature

### 3.1. add_*/mod_* functions

add functions are meant to add elements into the data, this could range from columns in a matrix to adding plots to the results section
mod functions on the other hand are meant to modify the existing elements in the data not that mod_mutate does add new columns to the data but we see it more as a modification performed on existing data and thus classified it as a mod function. 

```{r, eval=F}

setGeneric("[add|mod]_new_name", function(object, ...) standardGeneric("[add|mod]_new_name"))

setMethod("[add|mod]_new_name", "metime_analyser", function(object, ...) {
            # conditions to check if the arguments if an argument fails return the object as is
            # An example is given here
            if(!all(which_data %in% names(object@list_of_data))) {
                warning("dataset not found in the object. Exiting without making any changes")
                return(object)
            }

            # Add Logic to add elements to the object

            # finally update the functions_applied section of the results
            out <- object
            out <- add_function_info(object=out, function_name="[add|mod]_new_name", params=list(...))
            return(out)
    })

```

In the above example object is an instance of class metime_analyser. 

### 3.2. calc_*/meta_* function

calc functions are meant to perform methods that were discussed above. Similarly, meta functions are meant to perform meta-analysis on these calculations. Note that meta functions can only be applied onto the objects that contains the similar type of results. 

```{r, eval=F}

setGeneric("calc_new_name", function(object, which_data, cols_for_meta=NULL, name="calc_new_name_1", stratifications=NULL, ...) standardGeneric("calc_new_name"))

setMethod("calc_new_name", "metime_analyser", function(object, which_data, cols_for_meta=NULL, name="calc_new_name_1", stratifications=NULL, ...) {
            # conditions to check if the arguments if an argument fails return the object as is
            # An example is given here, of course there might be many more based on arguments.
            if(!all(which_data %in% names(object@list_of_data))) {
                warning("dataset not found in the object. Exiting without making any changes")
                return(object)
            }

            # Updating the name of the result here. Make sure to check the regex based on the name of the function
            if(grep(name, names(object@results)) %>% length() >=1) {
                warning("name of the results was previously used, using a different name")
                index <- name %>% gsub(pattern="[a-z|A-Z]+_[a-z|A-Z]+_[a-z|A-Z]+_", replacement="") %>% as.numeric()
                index <- c(0:9)[grep(index, 0:9)+1]
                name <- name %>% gsub(pattern="_[0-9]", replacement=paste("_", index, sep=""))
            }

            # Performing stratification analysis before proceeding with calculation
            data_list <- get_stratified_data(object=object, which_data=which_data, stratifications=stratifications)
            data <- data_list[["data"]]
            row_data <- data_list[["row_data"]]

            # get metadata based on cols_for_meta. These columns will later be used as color or shape aesthetics while plotting the results. Depending on the calculation choose metadata for metabolites or metadata for samples using get_metadata_for_columns() and get_metadata_for_rows() respectively. An example with metabolites is shown here
            metadata <- get_metadata_for_columns(object=object, which_data=which_data, columns=cols_for_meta)

            # Add Logic for calculation on data and row_data generated above here

            # Update results and function information here
            out <- get_make_results(object=object, data = result_here # this is a list, 
                                metadata = metadata, 
                                calc_type = rep("name_of_calc", each=length(result_here)), 
                                calc_info = "More detailed information as a string here",
                                name=name) %>%
                    add_function_info(object=out, function_name="calc_new_name", 
                            params=list(which_data=which_data, ..., cols_for_meta=cols_for_meta, 
                            name=name, stratifications=stratifications))
            return(out)
    })



```

In order to perform a calculation specific to a single dataset but want to do it for multiple datasets at once then wrap the code chunk from stratification to add_function_info() part using a loop. Also, make sure that the column names mentioned in stratifications list is consistent across the row_data of the datasets. To see such an example look at calc_conservation_* functions.

```{r, eval=F}

setGeneric("meta_new_name", function(object, results_index=NULL, name="meta_new_name_1", ...) standardGeneric("meta_new_name"))
setMethod("meta_new_name")

```
